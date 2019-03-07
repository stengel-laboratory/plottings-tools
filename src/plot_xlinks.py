#!/usr/bin/env python3.6

import argparse
import pandas as pd
import link_library as ll
import link_library.xtract_library.plot_xt as plot_xt

desc = """Kai Kammer - 2018-09-17. 
Script to visualize crosslink log2ratio patterns as given by xTract
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input csv files separated by spaces")
parser.add_argument('-o', '--outname', action="store", dest="outname", default='',
                    help="Append this postfix to the plot file name.")
args = parser.parse_args()


xt_db = ll.xTractDB()


@ll.timeit
def get_matching_monos(df):
    print(df.groupby(xt_db.type_string)[xt_db.uxid_string].nunique())
    def find_associated_link(x, df_monos):
        def is_equal(a,b):  # given two links a and b check whether they link the same protein and position
            prots_a = (a[xt_db.prot_string].iloc[0])
            pos_a = (a[xt_db.pos_string].iloc[0])
            prot_pos_a = set(zip(prots_a, pos_a))  # put prot name and pos into tuple for easy set intersection
            prots_b = (b[xt_db.prot_string].iloc[0])
            pos_b = (b[xt_db.pos_string].iloc[0])
            prot_pos_b = set(zip(prots_b, pos_b))
            link_intersect = prot_pos_a & prot_pos_b  # get intersecting links
            if len(link_intersect) > 0:
                return True
            return False

        tmp = df_monos.groupby(xt_db.uxid_string).filter(lambda y: is_equal(x, y))
        if len(tmp) > 0:
            tmp[xt_db.link_group_string] = x[xt_db.link_group_string].iloc[0]
            # x[associated_link_string] = [tmp[uid_string].values]  # assignment is buggy and not needed anyway
            x = pd.concat([x,tmp], sort=True)
        return x

    df_xlinks = df[df[xt_db.type_string] == xt_db.type_xlink_string].copy()  # using a copy since I set values in the next line
    df_xlinks[xt_db.link_group_string] = range(len(df_xlinks))
    df_monos = df[df[xt_db.type_string] == xt_db.type_mono_string]
    df_new = df_xlinks.groupby(xt_db.uxid_string).apply(lambda x: find_associated_link(x, df_monos)).reset_index(drop=True)
    num_mono1 = df_new.groupby(xt_db.link_group_string).filter(lambda x: len(x[x[xt_db.type_string] == xt_db.type_mono_string]) >= 1 and len(x[x[xt_db.type_string] == xt_db.type_xlink_string]) == 1)[xt_db.link_group_string].nunique()
    num_mono2 = df_new.groupby(xt_db.link_group_string).filter(lambda x: len(x[x[xt_db.type_string] == xt_db.type_mono_string]) == 2 and len(x[x[xt_db.type_string] == xt_db.type_xlink_string]) == 1)[xt_db.link_group_string].nunique()
    num_total = df_new[xt_db.link_group_string].nunique()
    print("Link groups with 1 monolink: {0} ({1:.0%})".format(num_mono1, num_mono1/num_total))
    print("Link groups with 2 monolinks: {0} ({1:.0%})".format(num_mono2, num_mono2/num_total))
    return df_new


def plot_mono_vs_xlink_quant(df, mono_1_plotter, mono_2_plotter):
    df_1 = filter_link_groups(df, 1)
    df_2 = filter_link_groups(df, 2)
    mono_1_plotter.plot_mono_vs_xlink_quant(df_1)
    mono_2_plotter.plot_mono_vs_xlink_quant(df_2)


def plot_associated_mono_links(df, df_dist, mono_1_plotter, mono_2_plotter):
    df_1 = filter_link_groups(df, 1)
    df_2 = filter_link_groups(df, 2)
    mono_1_plotter.plot_associated_mono_links(df_1, df_dist)
    mono_2_plotter.plot_associated_mono_links(df_2, df_dist)

def plot_dist_vs_quant(df, df_dist):
    plotter = plot_xt.PlotMaster('plots_quant')
    plotter.plot_dist_vs_quant(df, df_dist)


def renumber_groups(x):
    if not hasattr(renumber_groups, "counter"):
        renumber_groups.counter = 0  # it doesn't exist yet, so initialize it
    x[xt_db.link_group_string] = renumber_groups.counter
    renumber_groups.counter += 1
    return x


def rename_groups(x):
    uid_list = []
    entry_link = x[x[xt_db.type_string] == xt_db.type_xlink_string]
    uid_list.append(entry_link[xt_db.uxid_string].iloc[0])
    for n in range(len(x[x[xt_db.type_string] == xt_db.type_mono_string])):
        uid_list.append(entry_link[xt_db.uxid_string].iloc[0])
    x[xt_db.link_group_string] = uid_list
    return x


def filter_imputed(df):
    # == means link was found in both experiments; stupid libreoffice saves == as Err:520
    df = df[(df['sign'] == '==') | (df['sign'] == 'Err:520') ]
    return df


def filter_link_groups(df, no_monos):
    df = df.groupby(xt_db.link_group_string).filter(lambda x: len(x[x[xt_db.type_string] == xt_db.type_xlink_string]) == 1)
    df = df.groupby(xt_db.link_group_string).filter(lambda x: len(x[x[xt_db.type_string] == xt_db.type_mono_string]) >= no_monos)
    df = df.reset_index(drop=True)
    df = df.groupby(xt_db.link_group_string).apply(rename_groups)
    return df


def main():
    df_xtract = None
    df_dist = None
    for inp in args.input:
        if ".xls" in inp:
            # the xls files written by xtract are buggy and can only be read this way
            df = pd.read_csv(inp, delim_whitespace=True)
        else:
            df = pd.read_csv(inp, engine='python')
        if xt_db.uxid_string in df and not xt_db.dist_string in df:
            df_xtract = df
        elif xt_db.dist_string in df:
            df_dist = df
        else:
            print("WARNING: No compatible input found for {0}".format(inp))
    if df_dist is None:
        mono_1_plotter = plot_xt.PlotMaster("plots_quant/mono_1" + args.outname)
        mono_2_plotter = plot_xt.PlotMaster("plots_quant/mono_2" + args.outname)
    else:
        mono_1_plotter = plot_xt.PlotMaster("plots_quant/mono_1_dist" + args.outname)
        mono_2_plotter = plot_xt.PlotMaster("plots_quant/mono_2_dist" + args.outname)
    df_tmp = df_xtract[xt_db.uxid_string].str.split(':').apply(ll.get_prot_name_and_link_pos)
    # direct assignment does not work as it takes just the column headers as values
    df_xtract[xt_db.pos_string], df_xtract[xt_db.prot_string] = df_tmp[0], df_tmp[1]
    df_xtract = filter_imputed(df_xtract)
    df_xtract = get_matching_monos(df_xtract)
    df_xtract = df_xtract.sort_values([xt_db.link_group_string, xt_db.type_string])
    df_xtract.to_csv("link_groups_out.csv")
    plot_mono_vs_xlink_quant(df_xtract, mono_1_plotter, mono_2_plotter)
    plot_associated_mono_links(df_xtract, df_dist, mono_1_plotter, mono_2_plotter)
    if df_dist is not None:
        plot_dist_vs_quant(df_xtract, df_dist)


if __name__ == "__main__":
    main()
