#!/usr/bin/env python3.6

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker

desc = """Kai Kammer - 2018-09-17. 
Script to visualize crosslink patterns
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', action="store", dest="input", default=None, type=str, nargs='+',
                    help="List of input csv files separated by spaces")
parser.add_argument('-o', '--output', action="store", dest="output", default='plot.png',
                    help="Name for the output figure")
args = parser.parse_args()

uxid_string = 'uxID'
uid_string = 'uID'
pos1_string = 'AbsPos1'
pos2_string = 'AbsPos2'
score_string = 'ld-Score'
log2_string = 'log2ratio'
times_found_xq_string = 'times_found_xq'
fdr_string = 'FDR'
type_string = 'type'
type_xlink_string = 'xlink'
type_mono_string = 'monolink'
signi_string = 'significance'
mono_to_xlink_string = 'mono_to_xlink'
xlink_to_mono_string = 'xlink_to_mono'
link_group_string = 'link_group'
found_partners_string = 'partners_found'


def plot_links(df_xtract, dfs_xquest):
    # get dict like columns->values
    dict_xtract = df_xtract.to_dict(orient='list')
    uid_count_dict = {k: 0 for k in dict_xtract[uid_string]}
    for df_xquest in dfs_xquest:
        dict_xquest = df_xquest.to_dict(orient='list')
        for uid in dict_xtract[uid_string]:
            if uid in dict_xquest[uxid_string]:
                uid_count_dict[uid] += 1

    df_count = pd.DataFrame.from_dict(
        {uid_string: list(uid_count_dict.keys()), times_found_xq_string: list(uid_count_dict.values())},
        orient='columns')
    df = pd.merge(df_xtract, df_count, on=[uid_string], how='outer')
    sp = sns.scatterplot(x=uid_string, y=log2_string, hue=times_found_xq_string, size=type_string, data=df,
                         palette="Set2")
    sp.hlines(y=1, xmin=0, xmax=len(uid_count_dict), label="log2=1", colors=['purple'])
    sp.hlines(y=-1, xmin=0, xmax=len(uid_count_dict), label="log2=1", colors=['purple'])
    plt.tight_layout()
    # sns.barplot(x=uid_string,y=log2_string,data=df)
    # df.plot(x=uid_string,y=log2_string,kind='bar')
    # plt.xticks(rotation=30)
    # sp.set_xticklabels(sp.get_xticklabels(), rotation=45, ha='right')



def plot_associated_mono_links(df_xtract):
    df_xtract = get_matching_monolinks(df_xtract)
    df_xlink = df_xtract.loc[df_xtract[type_string] == type_xlink_string]
    df_xlink = df_xlink.loc[df_xlink[found_partners_string] != 0]
    # print(df_xlink)
    # check to see if lists are empty
    #df_xtract.loc[df_xtract[link_group_string].map(lambda d: len(d)) == 0] = 0
    # df_xtract.loc[df_xtract[found_partners_string].map(lambda d: len(d)) == 0] = 0
    df_xtract = df_xtract.loc[df_xtract[found_partners_string] > 0]
    df_xtract = df_xtract.loc[df_xtract[fdr_string] <= 0.05]
    ax = sns.scatterplot(x=link_group_string, y=log2_string, style=type_string, hue=type_string, data=df_xtract,
                        palette="Set1", s=150)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    b_alt = True
    for x in df_xtract[link_group_string]:
        if b_alt:
            c = 'lightcoral'
            b_alt = False
        else:
            c = 'skyblue'
            b_alt = True
        ax.vlines(x=x,ymin=min(df_xtract[log2_string]),ymax=max(df_xtract[log2_string]),linestyles=':',colors=c)
    plt.tight_layout()


def get_matching_monolinks(df_xtract: pd.DataFrame):
    df_new = pd.DataFrame(df_xtract)
    df_new = df_new.assign(mono_to_xlink=pd.Series([[] for l in range(len(df_new))]),
                           xlink_to_mono=pd.Series([[] for l in range(len(df_new))]),
                           link_group=pd.Series([] for l in range(len(df_new))),
                           partners_found=pd.Series(0 for l in range(len(df_new))))
    df_mono = df_xtract.loc[df_xtract[type_string] == type_mono_string]
    df_xlink = df_xtract.loc[df_xtract[type_string] == type_xlink_string]
    i_lnk_grp = 0
    # iterrows is painfully slow; should use something else: https://stackoverflow.com/questions/24870953/does-iterrows-have-performance-issues
    # now using converted dict
    dict_xlink = df_xlink.to_dict(orient='index')
    dict_mono = df_mono.to_dict(orient='index')
    for entry_xlink in dict_xlink.values():
        i_lnk_grp += 1
        for entry_mono in dict_mono.values():
            if entry_mono[uid_string] in entry_xlink[uid_string]:
                df_new.loc[df_new[uid_string] == entry_mono[uid_string], mono_to_xlink_string] += [entry_xlink[uid_string]]
                df_new.loc[df_new[uid_string] == entry_xlink[uid_string], xlink_to_mono_string] += [entry_mono[uid_string]]
                df_new.loc[df_new[uid_string] == entry_mono[uid_string], link_group_string] += [i_lnk_grp]
                df_new.loc[df_new[uid_string] == entry_xlink[uid_string], link_group_string] = i_lnk_grp
                df_new.loc[df_new[uid_string] == entry_mono[uid_string], found_partners_string] += 1
                df_new.loc[df_new[uid_string] == entry_xlink[uid_string], found_partners_string] += 1
            else:
                df_new.loc[df_new[uid_string] == entry_xlink[uid_string], link_group_string] = i_lnk_grp
    df_mono = df_new.loc[df_new[type_string] == type_mono_string]
    df_xlink = df_new.loc[df_new[type_string] == type_xlink_string]
    df_new = pd.DataFrame()
    for i_mono, entry_mono in df_mono.iterrows():
        partners_found = entry_mono[found_partners_string]
        if partners_found > 0:
            for group_id in entry_mono[link_group_string]:
                new_entry_mono = pd.Series(entry_mono)
                new_entry_mono.loc[link_group_string] = group_id
                df_new = df_new.append([new_entry_mono])
            # print(entry_mono[found_partners_string])
    # print(df_new.loc[df_xtract[type_string] == type_mono_string])
    df_new = df_new.append(df_xlink)
    # print(df_new.loc[df_xtract[type_string] == type_xlink_string])
    return df_new


def main():
    dfs_xquest = []
    df_xtract = None
    for inp in args.input:
        if ".xls" in inp:
            # the xls files written by xtract are buggy and can only be read this way
            df = pd.read_csv(inp, delim_whitespace=True)
        else:
            df = pd.read_csv(inp, engine='python')
        if uid_string in df:
            df_xtract = df
        elif uxid_string in df:
            dfs_xquest.append(df)
        else:
            print("WARNING: No compatible input found for {0}".format(inp))
    if dfs_xquest:
        plot_links(df_xtract, dfs_xquest)
    else:
        plot_associated_mono_links(df_xtract)
    plt.savefig(args.output)
    plt.show()

if __name__ == "__main__":
    main()
