#! /usr/bin/env python
import polars as pl
import polars
import sys
import argparse
from collections import defaultdict

RANK_IDX = {
        "superkingdom": 0,
        "kingdom": 0,
        "phylum": 1,
        "class": 2,
        "order": 3,
        "family": 4,
        "genus": 5,
        "species": 6
    }


def load_gather_lineages(csv_path, *, rank="class"):
    df = polars.read_csv(csv_path)

    lin_idx = RANK_IDX[rank]

    def split_lineages(row):
        val = row["lineage"]
        if val is None:
            return "UNASSIGNED"

        return ";".join(val.split(';')[:lin_idx+1])

    df = df.with_columns([
        pl.struct(["lineage"])
        .map_elements(split_lineages, return_dtype=str)
        .alias("sublineage") ])
    return df


def elide(s, length=70):
    if len(s) > length:
        return s[:length-3] + '...'
    return s
    


def main():
    p = argparse.ArgumentParser()
    p.add_argument('csvs', nargs='+', help='output of tax annotate on gather')
    p.add_argument('-f', '--force', help='continue past loading errors',
                   action='store_true')
    p.add_argument('-r', '--rank', help='aggregate at this rank',
                   default="class")
    p.add_argument('-m', '--min-lineages', help='only display records with lineages >= this number; default=2',
                   type=int, default=2)
    p.add_argument('--match-cutoff', type=float, default=0.01,
                   help='only display matches above this cutoff')
    args = p.parse_args()

    total_count = 0
    multi_count = 0

    for csv in args.csvs:
        try:
            full_df = load_gather_lineages(csv, rank=args.rank)
        except:
            if args.force:
                continue
            raise

        # split dataframe on query_name
        df_by_query = full_df.partition_by("query_name", as_dict=True)

        # for each query, ask: does it belong to multiple lineages?
        for (query_name,), df in df_by_query.items():
            df2 = df.select(["match_name", "lineage", "sublineage", "f_unique_weighted"])

            lin_sum = defaultdict(float)
            full_lineages = defaultdict(list)
            for match_name, lin, sublin, f_unique_weighted in df2.iter_rows():
                if lin is None:
                    # @CTB error out or something...? lineage unassigned by
                    # sourmash tax.
                    # assert 0, (csv, query_name, match_name)
                    lin = "UNASSIGNED"
                lin_sum[sublin] += f_unique_weighted
                full_lineages[sublin].append(lin)

            n_hidden = 0
            if len(lin_sum) >= args.min_lineages:
                multi_count += 1
                total = 0.0
                print(elide(query_name))
                for sublin, frac in sorted(lin_sum.items(),
                                           key=lambda x: -x[1]):

                    # how many lineages w/in this sublin? if 1, display full.
                    n_lin = len(full_lineages[sublin])
                    if n_lin == 1:
                        sublin = full_lineages[sublin][0]
                    else:
                        sublin = f"{sublin} ({n_lin} species)"
                    total += frac
                    if frac >= args.match_cutoff:
                        print(f'\t{frac * 100:.1f}% {elide(sublin, 65)}')
                    else:
                        n_hidden += 1
                if total != 1.0:
                    if n_hidden:
                        print(f'\t** {n_hidden} matches below threshold; {(1 - total)*100:.1f}% unknown')
                    else:
                        print(f'\t** {(1 - total)*100:.1f}% unknown')

        total_count += len(df_by_query)


if __name__ == '__main__':
    main()
