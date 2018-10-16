#!/usr/bin/env python3.6

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# TODO: uid level is fine; however for doing the violations on uxid level they would have to be calculated before sum()
# TODO: remove regular container support; only use detials container and regular containers as a control

desc = """Kai Kammer - 2018-09-17. 
Script to plot xTract bag container ms1 areas
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--input', action="store", dest="input", default=None, type=str, nargs='+',
                    help="List of input csv files separated by spaces")
parser.add_argument('-l', '--level_ms1', action="store", dest="level", default='uID',
                    help="Level on which the ms1 intensities are summed. Either uID (peptide)"
                         "or uxID (protein).")
parser.add_argument('-t', '--title', action="store", dest="title", default='',
                    help="Optional title for the plot.")
parser.add_argument('-f', '--filter', action="store", dest="filter", default="",
                    help="Optionally specify a link type to filter for. Possible values: monolink, xlink")
parser.add_argument('-p', '--plot_type', action="store", dest="plot", default='scatter',
                    help="Type of plot. Possible values: scatter, bar, lh (light heavy), rep (replicates), rep_bar")
parser.add_argument('-o', '--output', action="store", dest="output", default='plot_plot_type.png',
                    help="Name for the output figure")
args = parser.parse_args()

setting_filter = ""  # filter is set in __main__ if given

col_level = ""  # sets actual level for run in __main__ (i.e. uid or uxid)

col_uid = "uid"
col_uxid = "uxid"
col_exp = "exp_name"
col_exp_original = 'exp_name_original'
col_origin = 'origin'
col_link_type = 'link_type'
col_weight_type = 'weight_type'
col_bag_type = 'bag_container_type'
col_area_sum_total = 'ms1_area_sum'
col_area_sum_norm_total = col_area_sum_total + '_norm'
col_area_sum_light = 'ms1_area_sum_light'
col_area_sum_heavy = 'ms1_area_sum_heavy'
col_area_bio_repl = 'ms1_area_sum_bio_rep_'
col_var = 'ms1_area_variance'
col_std = 'ms1_area_std'
col_index = 'index'
col_log2ratio = 'log2ratio'
col_lh_log2ratio = 'light_heavy_log2ratio'
col_bio_rep = 'exp_bio_rep'
row_monolink_string = 'monolink'
row_xlink_string = 'xlink'
row_light_string = 'light'
row_heavy_string = 'heavy'
row_regular_string = 'regular'
row_details_string = 'details'

class BagContainer(object):
    def __init__(self, mode):
        if mode == row_regular_string:
            self.cont_type = row_regular_string
            self.uxid_string = 'a_uxID'
            self.uid_string = 'a_uID'
            self.exp_string = 'a_experimentname'
            self.vio_string = 'c_violations'
            self.sum_light_string = 'b_light_msum_area_sum_isotopes'
            self.sum_heavy_string = 'b_heavy_msum_area_sum_isotopes'

        else:
            self.cont_type = row_details_string
            self.uxid_string = 'b_peptide_uxID'
            self.uid_string = 'b_peptide_uID'  # use sequence instead of 'b_peptide_uID' as the latter contains too much information
            self.seq_string = 'b_peptide_seq'
            self.exp_string = 'd_exp_name'
            self.fraction_string = 'd_exp_fraction'
            self.repl_bio_string = 'd_exp_biol_rep'
            self.repl_tech_string = 'd_exp_tech_rep'
            self.vio_string = 'a_bag_container_violations'
            self.sum_string = 'c_pg_area_sum_isotopes'
            self.valid_string = "b_peptide_var_valid"
            self.type_string = "b_peptide_type"
        if col_level == col_uid:
            self.level = self.uid_string
        else:
            self.level = self.uxid_string


def plot_bag_container_bar(df):
    # filter zero ms1 intensities
    df = df.loc[df[col_area_sum_total] > 0]
    # filter ids not found in at least two experiments
    df = df.groupby(col_level).filter(lambda x: len(x) > 1)
    # filter by log2ratio
    # df = df[(df[col_log2ratio] > 3) | (df[col_log2ratio] < -3)]
    df = df.reset_index()
    # following line can be used to only plot partial data
    # df = df.loc[df[index_string] < len(df.index)]
    ax = sns.barplot(x=col_level, y=col_area_sum_total, hue=col_exp, data=df)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right', size=9)
    ax.set(title="{0} ({1} level)".format(args.title, args.level), yscale='log')
    save_n_show_fig(args.plot)

def plot_bag_container_scatter(df):
    exp_list = sorted(list(set(df[col_exp])))
    if len(exp_list) > 2:
        print("More than two experiments found. Please select which ones to plot.")
        print("{0}".format({no:exp for no,exp in enumerate(exp_list)}))
        exp1 = int(input("Please select first experiment: "))
        exp2 = int(input("Please select second experiment: "))
        exp1, exp2 = exp_list[exp1], exp_list[exp2]
    elif len(exp_list) == 2:
        exp1, exp2 =exp_list[0], exp_list[1]
    else:
        print("ERROR: Too few experiments: {0}".format(exp_list))
        exit(1)
    #df = df.loc[df[col_area_sum_norm_total] > 0]  # filter zero intensities

    df_x = df.loc[df[col_exp] == exp1]
    df_y = df.loc[df[col_exp] == exp2]
    df = pd.merge(df_x, df_y, on=[col_level, col_link_type], how='inner')  # inner: only merge intersection of keys
    df = df.dropna()
    df = df.reset_index()
    # df = df.loc[df[index_string] < len(df.index)]

    # note that regplot (underlying lmplot) will automatically remove zero values when using log scale
    fg = sns.lmplot(x=col_area_sum_total + '_x', y=col_area_sum_total + '_y', hue=col_link_type, data=df,
                    fit_reg=False, robust=False, ci=None)
    min_x = df[col_area_sum_total + '_x'].min()
    min_y = df[col_area_sum_total + '_y'].min()
    min_min = min(min_x, min_y)
    # using same minimum value for x and y and offset it by half to not cut off values at the limits
    min_min -= min_min/2
    fg.set(xlabel="{0} ({1})".format(col_area_sum_total,
                                     df[col_exp+'_x'][0]),  #KK_26S_merged_final.analyzer.quant
           ylabel="{0} ({1})".format(col_area_sum_total,
                                     df[col_exp+'_y'][0]),
           xscale='log', yscale='log', title="{0} ({1} level)".format(args.title, args.level),
           xlim=min_min,ylim=min_min)
    # draw horizontal line for all possible plots
    for row in fg.axes:
        for ax in row:
            ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
    save_n_show_fig(args.plot)


def plot_bag_container_light_heavy_scatter(df):
    # df_new = df_new[df_new[col_link_type] == row_xlink_string]
    fg = sns.lmplot(x=col_area_sum_light, y=col_area_sum_heavy, hue=col_link_type,
                    col=col_exp_original, row=col_origin, data=df, fit_reg=False, sharex=True, sharey=True, robust=True, ci=None, legend_out=False, )
    df_new = df[(df[col_area_sum_heavy] > 0) & (df[col_area_sum_light] > 0)]
    min_val = np.min(df_new[[col_area_sum_light, col_area_sum_heavy]].min())
    min_val -= min_val/2
    max_val = np.max(df_new[[col_area_sum_light, col_area_sum_heavy]].max())
    max_val += max_val/2
    # note that not setting x,ylim to auto (the default) leads to strange scaling bugs with a log scale
    # therefore using the same limits for all subplots; also makes comparisons easier
    fg.set(xscale='log', yscale='log',xlim=(min_val,max_val), ylim=(min_val,max_val))
    # draw horizontal line for all possible plots
    for row in fg.axes:
        for ax in row:
            ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
    save_n_show_fig(args.plot)


def plot_bag_container_bio_rep_scatter(df):
    bio_rep_list = [1,2,3]#sorted(list(set(df[col_area_bio_repl])))
    bio_rep_cols_list = [col_area_bio_repl+str(x) for x in bio_rep_list]
    # df_new = df_new[df_new[col_link_type] == row_xlink_string]
    for n_outer, bio_rep_outer in enumerate(bio_rep_cols_list):
        for n_inner, bio_rep_inner in enumerate(bio_rep_cols_list):
            if n_inner > n_outer:
                fg = sns.lmplot(x=bio_rep_outer, y=bio_rep_inner, #hue=col_link_type,
                                row=col_exp_original, data=df, fit_reg=False, sharex=True, sharey=True,
                                ci=None, legend_out=False, hue=col_link_type)
                df_new = df[df[col_area_sum_total] > 0]
                min_val = np.min(df_new[[col_area_sum_total]].min())
                min_val -= min_val/2
                max_val = np.max(df_new[[col_area_sum_total]].max())
                max_val += max_val/2
                # note that not setting x,ylim to auto (the default) leads to strange scaling bugs with a log scale
                # therefore using the same limits for all subplots; also makes comparisons easier
                fg.set(xscale='log', yscale='log',xlim=(min_val,max_val), ylim=(min_val,max_val))
                # draw horizontal line for all possible plots
                for row in fg.axes:
                    for ax in row:
                        ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
                save_n_show_fig(args.plot + "_" + str(n_outer+1) + "_vs_" + str(n_inner+1))


def plot_bag_container_bio_rep_bar(df):
    # filter zero ms1 intensities
    # df = df.loc[df[col_area_sum_total] > 0]
    # filter ids not found in at least two experiments
    # df = df.groupby(col_level).filter(lambda x: len(x) > 1)
    # filter by log2ratio
    # df = df[(df[col_log2ratio] > 3) | (df[col_log2ratio] < -3)]
    # df = df.reset_index()
    # following line can be used to only plot partial data
    df = df.reset_index()
    df[col_index] = df[col_index].astype(int)
    df = df.loc[df[col_index] < len(df.index)/4]
    fg = sns.catplot(kind="bar", x=col_level, y=col_area_sum_total, hue=col_area_bio_repl, data=df, row=col_exp_original, ci=None)
    fg.set(yscale='log')
    for row in fg.axes:
        for ax in row:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right', size=9)
    save_n_show_fig(args.plot)


def rename_columns(df, bag_cont):
    if bag_cont.cont_type == row_details_string:
        df = df.rename(index=str, columns={bag_cont.exp_string: col_exp_original,bag_cont.level: col_level,
                                           bag_cont.sum_string: col_area_sum_total,
                                           bag_cont.repl_bio_string: col_bio_rep,})
    return df


# function to get ms1 intensities grouped by experiment; works with bag_container.stats and .details
# returns new dataframe
def get_sum_ms1_intensities_df(df, bag_cont):
    # we split up the input by experiment and compute all ms1 intensities separately
    exp_list = list(set(df[bag_cont.exp_string]))
    # this list stores the dataframes by experiment
    df_exp_list = []

    # let's remove the weight from uid string; so we have the same uids for both container types
    if bag_cont.cont_type == row_details_string and col_level == col_uid:
        # using regex to match either :heavy or :light and all the following string (.*)
        df[bag_cont.uid_string] = df[bag_cont.uid_string].str.replace(":({0}|{1}).*"
                                                                      .format(row_light_string, row_heavy_string),"")
    # all ids (either uxid or uid) are put into one set and used for constructing the initial dataframe
    # this also allows for detection of ids missing from one of the two experiments
    id_list = list(set(df[bag_cont.level]))
    # iterating by experiment
    for exp in exp_list:
        # since we use the experiment to group the plotting, the exp_new name will also include the original file name
        exp_new = "{0} ({1}): {2}".format(exp, bag_cont.cont_type, df.name[:df.name.rfind('_')])
        # creating the results dataframe
        df_res = pd.DataFrame()
        kwargs = {col_level: id_list,
                  col_exp: pd.Series([exp_new for l in range(len(id_list))]),
                  col_exp_original: pd.Series([exp for l in range(len(id_list))]),
                  col_bag_type: pd.Series([bag_cont.cont_type for l in range(len(id_list))]),
                  col_origin: pd.Series([df.name for l in range(len(id_list))])}
        df_res = df_res.assign(**kwargs)

        # filtering the input dataframe by experiment name
        df_exp = df.loc[df[bag_cont.exp_string] == exp]
        # processing bag_container.stats
        if bag_cont.cont_type == row_regular_string:
            # first filtering violations found by xTract
            df_exp = df_exp.loc[df_exp[bag_cont.vio_string] == 0]
            # putting the link type (mono/xlink) into its own column
            df_exp[col_link_type] = df_exp[bag_cont.uid_string].str.split('::').str[2]
            # selecting the columns necessary to sum up ms1 intensities which are separated into heavy and light
            df_ms1_tot = df_exp[[bag_cont.level, bag_cont.sum_light_string, bag_cont.sum_heavy_string, col_link_type]]
            # renaming the columns before merging them into our results df
            df_ms1_tot = df_ms1_tot.rename(index=str, columns={bag_cont.level: col_level,
                                                               bag_cont.sum_light_string: col_area_sum_light,
                                                               bag_cont.sum_heavy_string: col_area_sum_heavy})
            # computing the total ms1 sum by computing the light+heavy sum along the x-axis (i.e. row-wise)
            df_ms1_tot[col_area_sum_total] = df_ms1_tot[[col_area_sum_light, col_area_sum_heavy]].sum(axis=1)
            # this will sum up the intensities for the same ids;
            # i.e. does nothing at uid level, sums up uids at uxid level
            df_ms1_tot = df_ms1_tot.groupby(
                [col_level, col_link_type], as_index=False)[col_area_sum_total, col_area_sum_light, col_area_sum_heavy].sum()
            # merging ms1 intensities with our results df
            # outer join: union of keys; inner: intersection
            df_res = pd.merge(df_res, df_ms1_tot, on=[col_level], how='inner')
        # processing bag_container.details
        else:
            # filtering these two means we get exactly the same results as from the regular bag container
            # removes violations (but no violations are calculated for monolinks)
            df_exp = df_exp.loc[df_exp[bag_cont.vio_string] == 0]
            df_exp = df_exp.loc[df_exp[bag_cont.valid_string] == 1]

            df_exp = df_exp.rename(index=str, columns={bag_cont.level: col_level})
            # create two separate columns for link type and weight (i.e. heavy or light)
            df_exp[col_link_type], df_exp[col_weight_type] = df_exp[bag_cont.type_string].str.split(':').str
            # summing up the total ms1 sum for the same id; again only really sums at uxid level
            df_ms1_tot = df_exp.groupby(
                [col_level, col_link_type])[bag_cont.sum_string].sum().reset_index(name=col_area_sum_total)
            df_res = pd.merge(df_res, df_ms1_tot, on=[col_level], how='inner')
            # the intensities for light and heavy ids have to calculated explicitly for details container
            # we group by weight and sum up the itensities
            df_ms1_lh = df_exp.groupby(
                [col_level, col_weight_type])[bag_cont.sum_string].sum().reset_index(name=col_area_sum_total)
            # we want the heavy and light intensities as separate columns and therefore pivot the dataframe here
            df_ms1_lh = pd.pivot_table(df_ms1_lh, values=col_area_sum_total, index=[col_level], columns=col_weight_type)
            df_ms1_lh = df_ms1_lh.reset_index()
            df_ms1_lh = df_ms1_lh.rename(
                index=str, columns={row_light_string: col_area_sum_light, row_heavy_string: col_area_sum_heavy})
            df_res = pd.merge(df_res, df_ms1_lh, on=[col_level], how='inner')
            df_ms1_bio_rep = df_exp.groupby(
                [col_level, bag_cont.repl_bio_string])[bag_cont.sum_string].sum().reset_index(name=col_area_sum_total)
            df_ms1_bio_rep = pd.pivot_table(df_ms1_bio_rep, values=col_area_sum_total, index=[col_level], columns=bag_cont.repl_bio_string)
            df_ms1_bio_rep = df_ms1_bio_rep.reset_index()
            # df_ms1_bio_rep = df_ms1_bio_rep.fillna(-1)
            rep_list = sorted(set(df_exp[bag_cont.repl_bio_string]))
            rep_name_dict = {x: col_area_bio_repl + str(x) for x in rep_list}
            df_ms1_bio_rep = df_ms1_bio_rep.rename(index=str, columns=rep_name_dict)
            df_res = pd.merge(df_res, df_ms1_bio_rep, on=[col_level], how='inner')
            df_res[col_var] = df_res[list(rep_name_dict.values())].var(axis=1)
            df_res[col_std] = df_res[list(rep_name_dict.values())].std(axis=1)
        # not sure if we should fill up na here with 0
        # if we do it would mean an id not detected in one experiment gets a 0 intensity
        # df_res = df_res.fillna(0)
        # optionally filter for link type; would be more efficient to do this earlier
        if setting_filter:
            df_res = df_res.loc[df_res[col_link_type] == setting_filter]
        # computes a normalized by maximum ms1 intensity; not used atm
        df_res[col_area_sum_norm_total] = df_res[col_area_sum_total] / df_res[col_area_sum_total].max()
        # computes the light/heavy log2ratio
        df_res[col_lh_log2ratio] = np.log2(df_res[col_area_sum_light] / df_res[col_area_sum_heavy])

        # df_res[col_area_sum_norm_total] = (df_res[col_area_sum_total] - df_res[col_area_sum_total].min()) / (df_res[col_area_sum_total].max() - df_res[col_area_sum_total].min())
        df_exp_list.append(df_res)
    df_final = pd.DataFrame()
    for dfl in df_exp_list:
        df_final = df_final.append(dfl)
    # computing the same violations as extract and removing them; right now only works on uid level
    df_final = get_violation_removed_df(df_final, bag_cont)

    # computing the log2ratio between the two experiments; hardcoded right now; at least reference should be a user setting
    df_log2 = pd.pivot_table(df_final, values=col_area_sum_total, index=[col_level], columns=col_exp).reset_index()
    df_log2 = df_log2.dropna()
    df_log2[col_log2ratio] = np.log2(df_log2.iloc[:,2]/df_log2.iloc[:,1])
    df_final = pd.merge(df_final, df_log2[[col_level,col_log2ratio]], on=[col_level], how='left')
    return df_final


# removes uids with violations; only checks for light/heavy log2ratio violations
# TODO: implement violation which requires all charge states of a peptide to be present in both experiments
# TODO: split up into violation assign (use separate columns) and filtering
# TODO: keep violation columns from xTract for comparison
def get_violation_removed_df(df, bag_cont):
    name = set(df[col_origin])
    print("Shape of {0} before filtering zero intensities: {1}.".format(name, df.shape))
    df = df[(df[col_area_sum_light] > 0) | (df[col_area_sum_heavy] > 0)]
    print("Shape of {0} before filtering lh log2ratio: {1}.".format(name, df.shape))
    df = df[(df[col_lh_log2ratio] < 1) & (df[col_lh_log2ratio] > -1)]
    print("Shape of {0} after filtering: {1}.".format(name, df.shape))
    # ms1 intensities have to be divided by tech_replicates*bio_replicates when coming from a details container
    # TODO: don't hardcode the number but get it from the input file
    if bag_cont.cont_type == row_details_string:
        df[[col_area_sum_total, col_area_sum_light, col_area_sum_heavy]] = \
            df[[col_area_sum_total, col_area_sum_light, col_area_sum_heavy]].apply(lambda x: x / 6)
    return df


def main():
    global col_level, setting_filter
    df_list = [] #  tuple (dataframe, bag_container)
    if args.level.lower() == col_uid:
        col_level = col_uid
    elif args.level.lower() == col_uxid:
        col_level = col_uxid
    else:
        print("ERROR: Improper ms1 level entered: {0}".format(args.level))
        exit(1)
    if args.filter.lower() == row_monolink_string:
        setting_filter = row_monolink_string
    elif args.filter.lower() == row_xlink_string:
        setting_filter = row_xlink_string
    elif args.filter.lower():
        print("ERROR: Improper link filter entered: {0}".format(args.filter))
        exit(1)
    con_details = BagContainer(row_details_string)
    con_reg = BagContainer(row_regular_string)
    for inp in args.input:
        if ".xls" in inp:
            # the xls files written by xtract are buggy and can only be read this way
            df = pd.read_csv(inp, engine='python', delimiter='\t', na_values=['-'])
        else:
            df = pd.read_csv(inp, engine='python')
        df.name = os.path.basename(inp)
        if con_reg.uid_string in df:
            df_list.append((df,con_reg))
        elif con_details.uid_string in df:
            df_list.append((df,con_details))
        else:
            print("WARNING: No compatible input found for {0}".format(args.input))
            exit(1)
    df_new = pd.DataFrame()
    for entry in df_list:
        df_new = df_new.append(get_sum_ms1_intensities_df(entry[0], entry[1]))
    df_new = df_new.sort_values([col_origin, col_exp_original, col_link_type])
    if args.plot == 'scatter':
        plot_bag_container_scatter(df_new)
    elif args.plot == 'bar':
        plot_bag_container_bar(df_new)
    elif args.plot == 'lh':
        plot_bag_container_light_heavy_scatter(df_new)
    elif args.plot == 'rep':
        plot_bag_container_bio_rep_scatter(df_new)
    elif args.plot == 'rep_bar':
        df_new = rename_columns(df_list[0][0], df_list[0][1])
        plot_bag_container_bio_rep_bar(df_new)
    else:
        print("WARNING: No compatible plot specified: {0}".format(args.input))
        exit(1)
    # calling tight layout as often as possible; only seems to improve the figure


def save_n_show_fig(out_name):
    plt.tight_layout()
    figure = plt.gcf()  # get current figure
    figure.set_size_inches(19,12)
    plt.tight_layout()
    plt.tight_layout()
    plt.tight_layout()
    plt.savefig("plot_" + out_name + ".png")
    plt.show()

if __name__ == "__main__":
    main()
