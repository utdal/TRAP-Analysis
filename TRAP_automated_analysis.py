"""
Targeted purification of polysomal mRNA (TRAP-Seq)
"""
# Author(s): Diana Tavares Ferreira(@dianatavf), Nikhil Nageshwar Inturi(@unikill066)
# Organization: The University of Texas at Dallas
# License: MIT License

# Imports
from scipy import stats
from trap import traplib  # Note: trap is a library used from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8019688/
import plotly.graph_objects as go
from log import logging_functionality
import plotly.express as px, warnings
from scipy.stats import percentileofscore
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import dendrogram, linkage
import pandas as pd, numpy as np, matplotlib as mpl, os
import matplotlib.pyplot as plt, seaborn as sns, itertools

pd.options.mode.chained_assignment = None
warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)


def plot_bar_chart(dataframe, group_list, prefix, counter, filepath):
    counter = counter + 1
    dataframe = dataframe[['Gene Name', 'Log2fc_{}_{}_{}'.format(group_list[1], group_list[0], prefix)]]
    dataframe.sort_values(by=['Log2fc_{}_{}_{}'.format(group_list[1], group_list[0], prefix)],
                          ascending=False, inplace=True)
    dataframe_head = dataframe.head(15) if all(
        dataframe.head(15)['Log2fc_{}_{}_{}'.format(group_list[1], group_list[0], prefix)]) >= 0 else dataframe[
        dataframe['Log2fc_{}_{}_{}'.format(group_list[1], group_list[0], prefix)] >= 0]
    dataframe_tail = dataframe.tail(15) if all(
        dataframe.head(15)['Log2fc_{}_{}_{}'.format(group_list[1], group_list[0], prefix)]) < 0 else dataframe[
        dataframe['Log2fc_{}_{}_{}'.format(group_list[1], group_list[0], prefix)] < 0]

    bar_plot_head = go.Bar(x=dataframe_head['Gene Name'].head(15),
                           y=dataframe_head['Log2fc_{}_{}_{}'.format(group_list[1], group_list[0], prefix)].head(15),
                           marker=dict(color=dataframe_head[
                               'Log2fc_{}_{}_{}'.format(group_list[1], group_list[0], prefix)].head(15),
                                       coloraxis="coloraxis"))
    bar_plot_tail = go.Bar(x=dataframe_tail['Gene Name'].tail(15),
                           y=dataframe_tail['Log2fc_{}_{}_{}'.format(group_list[1], group_list[0], prefix)].tail(15),
                           marker=dict(color=dataframe_tail[
                               'Log2fc_{}_{}_{}'.format(group_list[1], group_list[0], prefix)].tail(15),
                                       coloraxis="coloraxis"))

    bar_plot = go.Figure(data=[bar_plot_head, bar_plot_tail])
    bar_plot.update_xaxes(tickangle=-45)
    bar_plot.update_layout(title="Log2fc_{}_{}_{} vs Genes Plot".format(group_list[1], group_list[0], prefix),
                           plot_bgcolor='white',
                           xaxis=dict(showgrid=False, showticklabels=True, showline=True, linecolor='black'),
                           yaxis=dict(title='Log2fc value', showgrid=False, showticklabels=True, showline=True,
                                      linecolor='black'), showlegend=False, coloraxis=dict(colorscale='viridis'))

    temp_file_path = filepath + os.sep + '{}_bar_plot_expression_Log2fc_{}_{}_{}_'.format(
        counter, group_list[1], group_list[0], prefix)
    bar_plot.write_html(temp_file_path + 'tpm.html')
    bar_plot.write_image(temp_file_path + 'tpm.pdf')

    my_logger.info("Bar plot for differentially expressed genes is plotted here {}".format(
        filepath + os.sep + '{}_bar_plot_expression_Log2fc_{}_{}_{}_'.format(counter, group_list[1],
                                                                             group_list[0], prefix) + 'tpm.html'))
    return counter


def plot_kde(plots_list, plot_name_list, counter, filepath):
    for kde_plot, dend_name in zip(plots_list, plot_name_list):
        counter += 1
        traplib.kdensity(kde_plot)  # Call the def function from traplib
        kde_plot_temppath = f"{filepath}/{counter}_kde_{dend_name}_tpm.pdf"
        plt.savefig(kde_plot_temppath)
        my_logger.info("KDE plot saved to {}".format(kde_plot_temppath))
        plt.clf()  # Empty the plot object
    return counter


def plot_dendrogram(plots_list, plot_name_list, counter, filepath):
    for dend_plot, dend_name in zip(plots_list, plot_name_list):
        counter = counter + 1
        sns.set_style("ticks")
        samples = list(dend_plot.columns.values)
        D = np.transpose(dend_plot)
        Z = linkage(D, method='average', metric='correlation')
        plt.figure(figsize=(23, 12))
        ax = plt.subplot()
        plt.xlabel('Samples', fontsize=15)
        plt.ylabel('Distance', fontsize=15)
        plt.title("Dendrogram for {}".format(dend_name), fontsize=20)
        dendrogram(Z, labels=samples)
        plt.xticks(rotation=10, ha='right')
        dend_plot_temppath = filepath + os.sep + '{}'.format(counter) + '_dendrogram_' + dend_name + '_tpm.pdf'
        plt.savefig(dend_plot_temppath)
        my_logger.info("Dendrogram plot saved to {}".format(dend_plot_temppath))
        plt.clf()
    return counter


def plot_heatmaps(plots_list, plot_name_list, counter, filepath):
    for heatmp, dend_name in zip(plots_list, plot_name_list):
        counter += 1
        m = px.imshow(heatmp.corr(), color_continuous_scale="viridis", aspect=True)
        m.update_layout(title="Heatmap for {}".format(dend_name))
        heatmap_html_path = f"{filepath}/{counter}_heatmap_{dend_name}_tpm.html"
        heatmap_pdf_path = f"{filepath}/{counter}_heatmap_{dend_name}_tpm.pdf"
        my_logger.info("Heatmap is saved to {}".format(heatmap_html_path))
        m.write_html(heatmap_html_path)
        m.write_image(heatmap_pdf_path)
    return counter


def plot_genes_of_interest(dataframe, gene_list, sample2_group2_IP, sample1_group2_IN,
                           sample2_group1_IP, sample1_group1_IN, group_list, counter, filepath):
    genes_temp_df = dataframe[dataframe['Gene Name'].isin(gene_list)]
    enrichment_df2 = pd.DataFrame(genes_temp_df.loc[:, ['Gene Name']])
    enrichment_df2['LogFC_{}_IP_IN'.format(group_list[0])] = traplib.log2_fold_change(genes_temp_df[sample2_group2_IP],
                                                                                      genes_temp_df[sample1_group2_IN])
    enrichment_df2['LogFC_{}_IP_IN'.format(group_list[1])] = traplib.log2_fold_change(genes_temp_df[sample2_group1_IP],
                                                                                      genes_temp_df[sample1_group1_IN])
    fc_df = traplib.stack_fig(enrichment_df2, position=7, character='h', group_names=group_list)
    fc_df.sort_values('tpm', ascending=False, inplace=True)

    fig, ax = plt.subplots(figsize=(25, 15))
    sns.set_context("paper")
    sns.set(font_scale=2)
    sns.set_style("whitegrid")
    custom_palette = ("medium blue", "amber")
    fc = sns.barplot(x='Gene Name', y='tpm', hue='group', data=fc_df, ci=68, alpha=0.95,
                     capsize=0.08, palette=sns.xkcd_palette(custom_palette))
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xticks(rotation=15)
    fig.suptitle('Genes of Interest', fontsize=20)
    plt.xlabel('Gene Name', fontsize=15)
    plt.ylabel('TPM', fontsize=15)

    counter = counter + 1
    fig.savefig(filepath + os.sep + "{}_genes_of_interest_barplot_in_log_ratio.pdf".format(counter))
    my_logger.info("<---- Initial plot to check the orientation of genes of interest {} is saved here: {} ---->".format(
        gene_list, filepath + os.sep + "{}_genes_of_interest_barplot_in_log_ratio.pdf".format(counter)))
    plt.clf()
    return counter


def plot_scatter(plots_list, plot_name_list, n, counter, filepath):
    for heatmp, dend_name in zip(plots_list[:2], plot_name_list[:2]):
        scatter_plot_group1_data = heatmp.iloc[:, :n]
        scatter_plot_group2_data = heatmp.iloc[:, n:]
        plot_group1_sc_before = list(itertools.combinations(scatter_plot_group1_data.columns, 2))
        plot_group2_sc_before = list(itertools.combinations(scatter_plot_group2_data.columns, 2))

        for i in plot_group1_sc_before:
            counter += 1
            pear = stats.pearsonr(scatter_plot_group1_data[i[0]], scatter_plot_group1_data[i[1]])
            ax = sns.regplot(x=scatter_plot_group1_data[i[0]], y=scatter_plot_group1_data[i[1]],
                             data=scatter_plot_group1_data, scatter_kws={"s": 4},
                             line_kws={'label': "pearsonr = {0:.1f}".format(pear[0]), 'color': 'C1'})
            ax.legend()
            plot_temppath = f"{filepath}/{counter}_scatter_pearsonr_{dend_name}_tpm.pdf"
            my_logger.info("Scatter plot is saved to {}".format(plot_temppath))
            plt.savefig(plot_temppath)
            plt.clf()

        for i in plot_group2_sc_before:
            counter += 1
            pear = stats.pearsonr(scatter_plot_group2_data[i[0]], scatter_plot_group2_data[i[1]])
            ax = sns.regplot(x=scatter_plot_group2_data[i[0]], y=scatter_plot_group2_data[i[1]],
                             data=scatter_plot_group2_data, scatter_kws={"s": 4},
                             line_kws={'label': "pearsonr = {0:.1f}".format(pear[0]), 'color': 'C1'})
            ax.legend()
            plot_temppath = f"{filepath}/{counter}_scatter_pearsonr_{dend_name}_tpm.pdf"
            my_logger.info("Scatter plot is saved to {}".format(plot_temppath))
            plt.savefig(plot_temppath)
            plt.clf()
    return counter


def compute_percentile_threshold(dataframe, dataframe_values, percentiles_list, filepath, counter, prefix):
    listofdicts = list()
    for pctl in percentiles_list:
        tempdict = {'percentile': pctl}
        for column in dataframe_values:
            tempdict[column] = dataframe_values[column].quantile(pctl)
        listofdicts.append(tempdict)

    pctl_df = pd.DataFrame(listofdicts)
    counter += 1
    thresh_filepath = os.path.join(filepath, "{}_percentile_threshold_{}_tpm_before.csv".format(counter, prefix))
    pctl_df.to_csv(thresh_filepath, index=False)

    # tpm_percentile_df = pd.DataFrame()
    for col in dataframe_values:
        colvals = dataframe_values[col]
        percentile_of_val = [percentileofscore(colvals, val) / 100 for val in colvals]
        dataframe[col + '_percentile'] = percentile_of_val
    print("")
    print("Check the {} percentile dataframe before setting the {} threshold: {}".format(
        prefix, prefix, thresh_filepath))
    my_logger.info("{} percentile dataframe is saved to {} before setting the {} threshold".format(
        prefix, thresh_filepath, prefix))

    return dataframe, thresh_filepath, counter


def convert_to_tpm(dataframe):
    numeric_cols = dataframe.select_dtypes(include=[float, int]).columns
    for column in numeric_cols:
        total_counts = dataframe[column].sum()
        if total_counts != 1000000:
            dataframe[column] = (dataframe[column] / total_counts) * 1000000
    return dataframe


def trap_analysis(sample1_group1_IN: list, sample2_group1_IP: list,
                  sample1_group2_IN: list, sample2_group2_IP: list,
                  filepath: str, gene_list_id: list, plot_name_list: list, group_list: list):
    """
        This method built around the analysis of TRAP data, given the mentioned parameters ...
            - Percentile calculation and filtering based on percentiles
            - Log2FC Calculations
            - Quantile Normalisation
            - Differential expression filter
        Additionally, several plots are created during the run-time, and are stored to "Processed_files" directory
        located at the same filepath as the input file.

        Example Samples:
            sample1_group1_IN = ['102_GAD65_IN_SNI_60_tpm', '104_GAD65_IN_SNI_60_tpm', '106_GAD65_IN_SNI_60_tpm']
            sample2_group1_IP = ['101_GAD65_IP_SNI_60_tpm', '103_GAD65_IP_SNI_60_tpm', '105_GAD65_IP_SNI_60_tpm']
            sample1_group2_IN = ['202_GAD65_IN_SHAM_60_tpm', '204_GAD65_IN_SHAM_60_tpm', '206_GAD65_IN_SHAM_60_tpm']
            sample2_group2_IP = ['201_GAD65_IP_SHAM_60_tpm', '203_GAD65_IP_SHAM_60_tpm', '205_GAD65_IP_SHAM_60_tpm']

        Note:
            Please change this following in the code: plot_name_list = ['gad_60_IP_values', 'gad_60_IN_values',
            'gad_60_values'], search for the above and change the strings accordingly ...
    """
    n = len(sample1_group1_IN)
    counter, sf = 0, 0.01
    csv_filepath = filepath
    group_list = [i.upper() for i in group_list]
    filepath = os.path.dirname(filepath) + os.sep + "Processed_files"
    percentiles_list = [round(x, 2) for x in np.arange(0.1, 1, 0.05)]

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    # loading the data [Note: Only xlsx, csv and tab file formats are supported at the moment]
    print("\n\nLoading the data from {}\n".format(csv_filepath))
    if csv_filepath.endswith('.csv'):
        tpm_dataframe = pd.read_csv(csv_filepath)
    elif csv_filepath.endswith('.xlsx'):
        tpm_dataframe = pd.read_excel(csv_filepath)
    elif csv_filepath.endswith('.tab'):
        tpm_dataframe = pd.read_csv(csv_filepath, header=0, delimiter="\t")
    else:
        raise Exception("Unrecognized file type, only .csv, .xlsx and .tab file formats are supported.")
    print("\n\nThe initial data has {} rows and {} columns\n".format(tpm_dataframe.shape[0],
                                                                     tpm_dataframe.shape[1]))
    my_logger.info("<---- The initial data has {} rows and {} columns ---->".format(tpm_dataframe.shape[0],
                                                                                    tpm_dataframe.shape[1]))

    # Check and convert to TPM's
    tpm_dataframe = convert_to_tpm(tpm_dataframe)

    gene_name = ['Gene Name']
    sample1_columns = sample1_group1_IN + sample1_group2_IN
    sample2_columns = sample2_group1_IP + sample2_group2_IP
    sample1_and_2_columns = sample1_group1_IN + sample1_group2_IN + sample2_group1_IP + sample2_group2_IP
    tpm_g_IN_IP_dataframe = tpm_dataframe.loc[:, gene_name + sample1_and_2_columns]

    tpm_IN_dataframe_values = tpm_dataframe.loc[:, sample1_columns]
    tpm_IP_dataframe_values = tpm_dataframe.loc[:, sample2_columns]
    tpm_g_IN_IP_dataframe_values = tpm_dataframe.loc[:, sample1_and_2_columns]

    plots_list = [tpm_IP_dataframe_values, tpm_IN_dataframe_values, tpm_g_IN_IP_dataframe_values]
    my_logger.info("<---- Loading the data to respective dataframes ---->")

    print("\n<---------------- Plots before the TRAP Analysis --------------->\n")
    my_logger.info("<---- Plots before the TRAP Analysis ---->")

    # bar-plot for genes of interest
    counter = plot_genes_of_interest(tpm_dataframe, gene_list_id, sample2_group2_IP, sample1_group2_IN,
                                     sample2_group1_IP, sample1_group1_IN, group_name_list, counter, filepath)

    # dendrogram plot -- before any modifications
    sns.set(font_scale=1)
    counter = plot_dendrogram(plots_list, plot_name_list, counter, filepath)

    # kde plots -- before any modifications
    counter = plot_kde(plots_list, plot_name_list, counter, filepath)

    # heatmaps -- before any modifications
    counter = plot_heatmaps(plots_list, plot_name_list, counter, filepath)

    # scatter plot -- before any modifications [Reg Plot]
    counter = plot_scatter(plots_list, plot_name_list, n, counter, filepath)
    print("Plots generated before TRAP Analysis are saved to {} directory.".format(filepath))
    my_logger.info("Plots generated before TRAP Analysis are saved to {} directory.".format(filepath))

    print("\n<---------------- TRAP Analysis --------------->\n")
    my_logger.info("<---- TRAP Analysis ---->")

    # additional check(s)
    if not (tpm_g_IN_IP_dataframe is not None and
            tpm_g_IN_IP_dataframe_values is not None and
            tpm_IN_dataframe_values is not None and tpm_IP_dataframe_values is not None):
        my_logger.exception("Dataframe is not defined.")
        raise Exception("Dataframe is not defined.")

    # Percentile calculation for IN cols
    my_logger.info("------------------------------------------------------------")
    my_logger.info("<------------- Adding percentile cols for IN -------------->")
    my_logger.info("------------------------------------------------------------")

    counter = counter + 1
    traplib.kdensity(tpm_g_IN_IP_dataframe[sample1_columns])
    kde_plot_temppath = filepath + os.sep + '{}_kde_plot_before_IN_qn_norm.pdf'.format(counter)
    plt.savefig(kde_plot_temppath)
    plt.clf()

    while True:
        try:
            in_norm = input(
                "Check the IN(input) kernel distribution {} and input if normalization is needed (Y/N) : ".format(
                    kde_plot_temppath))
            if in_norm in ['Y', 'N']:
                break
            else:
                print("The IN normalization input value takes only Y/N.")
        except ValueError:
            print("The IN normalization input value takes only Y/N.")
            continue

    tpm_g_IN_IP_dataframe, thresh_filepath_in, counter = compute_percentile_threshold(
        tpm_g_IN_IP_dataframe, tpm_IN_dataframe_values, percentiles_list, filepath, counter, "IN")

    sample1_group1_percentile_columns = [i + '_percentile' for i in sample1_group1_IN]  # sample1_group1_IN_temp
    sample1_group2_percentile_columns = [i + '_percentile' for i in sample1_group2_IN]  # sample1_group2_IN_temp

    percentile_sample1_IN_dataframe = tpm_g_IN_IP_dataframe[
        sample1_group1_percentile_columns + sample1_group2_percentile_columns]  # df_gad_IN_60
    log2_in = pd.DataFrame(np.log2(percentile_sample1_IN_dataframe.clip(lower=0.01)))
    temp_cols = [col + '_log2' for col in log2_in.columns]
    log2_in.columns = temp_cols
    log2_group1_in = log2_in.iloc[:, 0:n]  # group1 IN
    log2_group2_in = log2_in.iloc[:, n:]  # group2 IN
    tpm_g_IN_IP_dataframe['ssmd_in_{}_{}'.format(group_list[1],
                                                 group_list[0])] = traplib.ssmd(log2_group1_in, log2_group2_in)

    if in_norm == 'Y':
        qn_in = traplib.quantileNormalize(tpm_g_IN_IP_dataframe[sample1_columns])
        sample1_qn_columns = [col + '_qn' for col in qn_in.columns]  # qn_col_names_in
        qn_in.columns = sample1_qn_columns
        tpm_g_IN_IP_dataframe = pd.concat([tpm_g_IN_IP_dataframe, qn_in], axis=1)  # keep everything together
    else:
        qn_in = tpm_g_IN_IP_dataframe[sample1_columns]
        sample1_qn_columns = [col + '_qn' for col in qn_in.columns]  # qn_col_names_in
        qn_in.columns = sample1_qn_columns
        tpm_g_IN_IP_dataframe = pd.concat([tpm_g_IN_IP_dataframe, qn_in], axis=1)  # keep everything together

    while True:
        try:
            threshold = float(input("Enter the threshold for IN values: "))
            if 0 <= threshold <= 1:
                break
            else:
                print("The IN threshold value should be in the range of 0-1 and not {}, please try again.".format(threshold))
        except ValueError:
            print("Try entering a float value (0.0 - 1.0)")
            continue

    my_logger.info("IN threshold set to {}".format(threshold))

    tpm_percentile_qn_IN_dataframe = tpm_g_IN_IP_dataframe.loc[
        (tpm_g_IN_IP_dataframe[sample1_group1_percentile_columns] >= threshold).all(axis=1) |
        (tpm_g_IN_IP_dataframe[sample1_group2_percentile_columns] >= threshold).all(axis=1)]

    input_threshold_df = tpm_percentile_qn_IN_dataframe  # IN filtered dataframe

    counter = counter+1
    tpm_percentile_qn_IN_dataframe.to_csv(filepath + os.sep + "{}_percentile_threshold_IN_percentile_qn_tpm_after.csv".format(counter), index=False)
    my_logger.info("IN percentile data frame is saved to {} after setting the threshold".format(
        filepath + os.sep + "{}_percentile_threshold_IN_percentile_qn_tpm_after.csv".format(counter)))

    counter = counter + 1
    traplib.kdensity(tpm_percentile_qn_IN_dataframe[sample1_qn_columns])
    kde_plot_temppath = filepath + os.sep + '{}_kde_plot'.format(counter) + '_qn_tpm.pdf'
    plt.savefig(kde_plot_temppath)
    plt.clf()
    my_logger.info("Check the density plot after calculating the first percentile for IN data: {}".format(
        kde_plot_temppath))

    # Percentile calculation for IP cols
    my_logger.info("------------------------------------------------------------")
    my_logger.info("<------------- Adding percentile cols for IP -------------->")
    my_logger.info("------------------------------------------------------------")
    print("")
    counter = counter + 1
    traplib.kdensity(tpm_g_IN_IP_dataframe[sample1_columns])
    kde_plot_temppath = filepath + os.sep + '{}_kde_plot_before_IP_qn_norm.pdf'.format(counter)
    plt.savefig(kde_plot_temppath)
    plt.clf()

    while True:
        try:
            ip_norm = input(
                "Check the IP kernel distribution {} and input if normalization is needed (Y/N) : ".format(
                    kde_plot_temppath))
            if ip_norm in ['Y', 'N']:
                break
            else:
                print("The IP normalization input value takes only Y/N.")
        except ValueError:
            print("The IP normalization input value takes only Y/N.")
            continue

    tpm_ip_df_values = tpm_percentile_qn_IN_dataframe.loc[:, sample2_columns]

    tpm_percentile_qn_IN_dataframe, thresh_filepath_ip, counter = compute_percentile_threshold(
        tpm_percentile_qn_IN_dataframe, tpm_ip_df_values, percentiles_list, filepath, counter, "IP")

    sample2_group1_percentile_columns = [i+'_percentile' for i in sample2_group1_IP]
    sample2_group2_percentile_columns = [i+'_percentile' for i in sample2_group2_IP]

    tpm_percentile_qn_IN_and_IP_dataframe = tpm_percentile_qn_IN_dataframe
    temp_log_df = tpm_percentile_qn_IN_and_IP_dataframe[sample2_group1_percentile_columns+sample2_group2_percentile_columns]
    log2_ip = pd.DataFrame(np.log2(temp_log_df.clip(lower=0.01)))
    temp_cols = [col + '_log2' for col in log2_ip.columns]
    log2_ip.columns = temp_cols
    log2_group1_ip = log2_ip.iloc[:, 0:n]
    log2_group2_ip = log2_ip.iloc[:, n:]
    tpm_percentile_qn_IN_and_IP_dataframe['ssmd_ip_{}_{}'.format(group_list[1],
                                                                 group_list[0])] = traplib.ssmd(log2_group1_ip, log2_group2_ip)

    if ip_norm == 'Y':
        qn_ip = traplib.quantileNormalize(tpm_percentile_qn_IN_and_IP_dataframe[sample2_columns])
        sample2_qn_columns = [col + '_qn' for col in qn_ip.columns]
        qn_ip.columns = sample2_qn_columns
        tpm_percentile_qn_IN_and_IP_dataframe = pd.concat([tpm_percentile_qn_IN_and_IP_dataframe, qn_ip], axis=1)
    else:
        qn_ip = tpm_percentile_qn_IN_and_IP_dataframe[sample2_columns]
        sample2_qn_columns = [col + '_qn' for col in qn_ip.columns]
        qn_ip.columns = sample2_qn_columns
        tpm_percentile_qn_IN_and_IP_dataframe = pd.concat([tpm_percentile_qn_IN_and_IP_dataframe, qn_ip], axis=1)

    while True:
        try:
            threshold = float(input("Enter the threshold for IP values: "))
            if 0 <= threshold <= 1:
                break
            else:
                print("The IP threshold value should be in the range of 0-1 and not {}, please try again.".format(threshold))
        except ValueError:
            print("Try entering a float value (0.0 - 1.0)")
            continue

    my_logger.info("IP threshold set to {}".format(threshold))

    tpm_percentile_qn_IN_and_IP_dataframe = tpm_percentile_qn_IN_and_IP_dataframe.loc[
        (tpm_percentile_qn_IN_and_IP_dataframe[sample2_group1_percentile_columns] >= threshold).all(axis=1) |
        (tpm_percentile_qn_IN_and_IP_dataframe[sample2_group2_percentile_columns] >= threshold).all(axis=1)]

    # plotting kde plot for QN
    counter = counter + 1
    traplib.kdensity(tpm_percentile_qn_IN_and_IP_dataframe[list(qn_ip.columns)])
    kde_plot_temppath = filepath + os.sep + '{}_kde_plot'.format(counter) + '_qn_tpm.pdf'
    plt.savefig(kde_plot_temppath)
    plt.clf()
    my_logger.info("\n Check the density plot after calculating qn, log and percentiles_list on the IP data: {}".format(
        kde_plot_temppath))

    # ##########################################################################################
    #                               STATISTICAL ANALYSIS BC
    # ##########################################################################################
    qn_ip_df_group1 = tpm_percentile_qn_IN_and_IP_dataframe.loc[:, sample2_qn_columns[:n]]
    qn_ip_df_group2 = tpm_percentile_qn_IN_and_IP_dataframe.loc[:, sample2_qn_columns[n:]]
    tpm_percentile_qn_IN_and_IP_dataframe['Log2fc_{}_{}_IP'.format(group_list[1],
                                                                   group_list[0])] = traplib.log2_fold_change(
        qn_ip_df_group1, qn_ip_df_group2, averaging='median')
    qn_in_df_group1 = tpm_percentile_qn_IN_and_IP_dataframe.loc[:, sample1_qn_columns[:n]]
    qn_in_df_group2 = tpm_percentile_qn_IN_and_IP_dataframe.loc[:, sample1_qn_columns[n:]]
    tpm_percentile_qn_IN_and_IP_dataframe['Log2fc_{}_{}_IN'.format(group_list[1],
                                                                   group_list[0])] = traplib.log2_fold_change(
        qn_in_df_group1, qn_in_df_group2, averaging='median')
    tpm_percentile_qn_IN_and_IP_dataframe['BC_{}_{}_IP'.format(group_list[0],
                                                               group_list[1])] = pd.DataFrame(traplib.bhatt_distance(
        qn_ip_df_group1, qn_ip_df_group2))
    tpm_percentile_qn_IN_and_IP_dataframe['BC_{}_{}_IN'.format(group_list[0],
                                                               group_list[1])] = pd.DataFrame(traplib.bhatt_distance(
        qn_in_df_group1, qn_in_df_group2))
    counter = counter + 1
    tpm_percentile_qn_IN_and_IP_dataframe.to_csv(
        filepath + os.sep + "{}_percentile_threshold_IP_percentile_qn_tpm_after.csv".format(counter), index=False)
    my_logger.info("IP percentile data frame is saved to {} after setting the threshold".format(
        filepath + os.sep + "{}_percentile_threshold_IP_percentile_qn_tpm_after.csv".format(counter)))

    # Differential Expression Filtering
    while True:
        try:
            print("")
            bc_var = float(input("Enter the filter value for BC: "))
            if 0.0 <= bc_var <= 10.0:
                break
            else:
                print("The BC value: {} entered is not in range. Input a correct value between [0.0, 10.0]".format(
                    bc_var))
        except ValueError:
            print("Try entering a float value (0.0 - 1.0)")
            continue

    while True:
        try:
            ssmd_var = float(input("Enter the filter value for SSMD: "))
            if 0.0 <= ssmd_var <= 10.0:
                break
            else:
                print("The SSMD value: {} entered is not in range. Input a correct value between [0.0, 10.0]".format(
                    ssmd_var))
        except ValueError:
            print("Try entering a float value (0.0 - 1.0)")
            continue

    while True:
        try:
            log2fc_var = float(input("Enter the filter value for Log2fc: "))
            if 0.0 <= log2fc_var <= 10.0:
                break
            else:
                print("The Log2fc value: {} entered is not in range. Input a correct value between [0.0, 10.0]".format(
                    log2fc_var))
        except ValueError:
            print("Try entering a float value (0.0 - 1.0)")
            continue

        # data-filter validation checks
        if not min(tpm_percentile_qn_IN_and_IP_dataframe[
                       'BC_{}_{}_IP'.format(group_list[0], group_list[1])]) <= bc_var <= max(
                tpm_percentile_qn_IN_and_IP_dataframe['BC_{}_{}_IP'.format(group_list[0], group_list[1])]):
            my_logger.exception(
                "The range of bc_var is {} and {}, and the value entered is {}, please enter a value in range".format(
                    min(tpm_percentile_qn_IN_and_IP_dataframe['BC_{}_{}_IP'].format(group_list[0], group_list[1])),
                    max(tpm_percentile_qn_IN_and_IP_dataframe['BC_{}_{}_IP'.format(group_list[0], group_list[1])]),
                    bc_var))
            raise Exception(
                "The range of bc_var is {} and {}, and the value entered is {}, please enter a value in range".format(
                    min(tpm_percentile_qn_IN_and_IP_dataframe['BC_{}_{}_IP'.format(group_list[0], group_list[1])]),
                    max(tpm_percentile_qn_IN_and_IP_dataframe['BC_{}_{}_IP'.format(group_list[0], group_list[1])]),
                    bc_var))
        if not min(tpm_percentile_qn_IN_and_IP_dataframe[
                       'ssmd_ip_{}_{}'.format(group_list[1], group_list[0])]) <= ssmd_var <= max(
                tpm_percentile_qn_IN_and_IP_dataframe['ssmd_ip_{}_{}'.format(group_list[1], group_list[0])]):
            my_logger.exception(
                "The range of ssmd_var is {} and {}, and the value entered is {}, please enter a value in range".format(
                    min(tpm_percentile_qn_IN_and_IP_dataframe['ssmd_ip_{}_{}'.format(group_list[1], group_list[0])]),
                    max(tpm_percentile_qn_IN_and_IP_dataframe['ssmd_ip_{}_{}'.format(group_list[1], group_list[0])]),
                    ssmd_var))
            raise Exception(
                "The range of ssmd_var is {} and {}, and the value entered is {}, please enter a value in range".format(
                    min(tpm_percentile_qn_IN_and_IP_dataframe['ssmd_ip_{}_{}'.format(group_list[1], group_list[0])]),
                    max(tpm_percentile_qn_IN_and_IP_dataframe['ssmd_ip_{}_{}'.format(group_list[1], group_list[0])]),
                    ssmd_var))
        if not min(tpm_percentile_qn_IN_and_IP_dataframe[
                       'Log2fc_{}_{}_IP'.format(group_list[1], group_list[0])]) <= np.log2(log2fc_var) <= max(
                tpm_percentile_qn_IN_and_IP_dataframe['Log2fc_{}_{}_IP'.format(group_list[1], group_list[0])]):
            my_logger.exception(
                "The range of bc_var is {} and {}, and the value entered is {}, please enter a value in range".format(
                    min(tpm_percentile_qn_IN_and_IP_dataframe['Log2fc_{}_{}_IP'.format(group_list[1], group_list[0])]),
                    max(tpm_percentile_qn_IN_and_IP_dataframe['Log2fc_{}_{}_IP'.format(group_list[1], group_list[0])]),
                    np.log2(log2fc_var)))
            raise Exception(
                "The range of bc_var is {} and {}, and the value entered is {}, please enter a value in range".format(
                    min(tpm_percentile_qn_IN_and_IP_dataframe['Log2fc_{}_{}_IP'.format(group_list[1], group_list[0])]),
                    max(tpm_percentile_qn_IN_and_IP_dataframe['Log2fc_{}_{}_IP'.format(group_list[1], group_list[0])]),
                    np.log2(log2fc_var)))

    my_logger.info("BC threshold set to {}".format(bc_var))
    my_logger.info("SSMD threshold set to {}".format(ssmd_var))
    my_logger.info("Log2FC threshold set to {}".format(log2fc_var))

    tpm_percentile_qn_IN_and_IP_dataframe_rev_de_filter = tpm_percentile_qn_IN_and_IP_dataframe
    tpm_percentile_qn_IN_and_IP_dataframe = tpm_percentile_qn_IN_and_IP_dataframe[
        (tpm_percentile_qn_IN_and_IP_dataframe['BC_{}_{}_IP'.format(group_list[0], group_list[1])].abs() >= bc_var) &
        (tpm_percentile_qn_IN_and_IP_dataframe[
             'ssmd_ip_{}_{}'.format(group_list[1], group_list[0])].abs() >= ssmd_var) &
        (tpm_percentile_qn_IN_and_IP_dataframe['Log2fc_{}_{}_IP'.format(group_list[1], group_list[0])].abs() >= np.log2(log2fc_var))]

    log2ipin_dataframe = tpm_percentile_qn_IN_and_IP_dataframe.loc[:, ['Gene Name']]
    log2ipin_dataframe['Log2FC_{}_IP_IN'.format(group_list[1])] = traplib.log2_fold_change(
        tpm_percentile_qn_IN_and_IP_dataframe[sample2_qn_columns[:n]],
        tpm_percentile_qn_IN_and_IP_dataframe[sample1_qn_columns[:n]])
    log2ipin_dataframe['Log2FC_{}_IP_IN'.format(group_list[0])] = traplib.log2_fold_change(
        tpm_percentile_qn_IN_and_IP_dataframe[sample2_qn_columns[n:]],
        tpm_percentile_qn_IN_and_IP_dataframe[sample1_qn_columns[n:]])

    # log 2 fc log-log
    epsilon = 0.01
    log2ipin_dataframe['Log2FC_{}_{}_IP_IN'.format(group_list[0], group_list[1])] = np.log2(
        tpm_percentile_qn_IN_and_IP_dataframe[sample2_qn_columns[:n]].median(axis=1) + epsilon /
        tpm_percentile_qn_IN_and_IP_dataframe[sample1_qn_columns[:n]].median(axis=1) + epsilon / (
        tpm_percentile_qn_IN_and_IP_dataframe[sample2_qn_columns[n:]].median(axis=1) + epsilon /
        tpm_percentile_qn_IN_and_IP_dataframe[sample1_qn_columns[n:]].median(axis=1) + epsilon))

    temp_log2ipin_dataframe = log2ipin_dataframe
    counter = counter + 1
    log2ipin_dataframe.to_csv(filepath + os.sep + '{}_log_IP_IN_data'.format(counter) + '.csv', index=False)

    log2ipin_dataframe = log2ipin_dataframe.sort_values('Log2FC_{}_IP_IN'.format(group_list[1]), ascending=False)
    if log2ipin_dataframe.shape[0] > 50:
        log2ipin_dataframe = log2ipin_dataframe.head(50)
        log2ipin_dataframe = log2ipin_dataframe[log2ipin_dataframe['Log2FC_{}_IP_IN'.format(group_list[1])] > 0]
    else:
        log2ipin_dataframe = log2ipin_dataframe[log2ipin_dataframe['Log2FC_{}_IP_IN'.format(group_list[1])] > 0]
    log_fig = px.histogram(log2ipin_dataframe, x="Gene Name", y="Log2FC_{}_IP_IN".format(group_list[1]))
    log_fig.update_layout(xaxis=dict(showgrid=False, showticklabels=True, showline=True, linecolor='black'),
                          yaxis=dict(showgrid=False, showticklabels=True, showline=True, linecolor='black'),
                          plot_bgcolor='white')
    counter = counter + 1
    log_fig.write_html(filepath+os.sep+"{}_IP_IN_barplot_gene_name.html".format(counter))
    log_fig.write_image(filepath + os.sep + "{}_IP_IN_barplot_gene_name.pdf".format(counter))
    my_logger.info("Barplot is saved to {}, only the top genes are displayed here".format(
        filepath+os.sep+"{}_IP_IN_barplot_gene_name.html".format(counter)))

    temp_log2ipin_dataframe = temp_log2ipin_dataframe.sort_values(
        'Log2FC_{}_{}_IP_IN'.format(group_list[0], group_list[1]), ascending=False)
    log2ipin_dataframe = temp_log2ipin_dataframe
    if log2ipin_dataframe.shape[0] > 50:
        log2ipin_dataframe = log2ipin_dataframe.head(50)
        log2ipin_dataframe = log2ipin_dataframe[log2ipin_dataframe['Log2FC_{}_{}_IP_IN'.format(group_list[0],
                                                                                               group_list[1])] > 0]
    else:
        log2ipin_dataframe = log2ipin_dataframe[log2ipin_dataframe['Log2FC_{}_{}_IP_IN'.format(group_list[0],
                                                                                               group_list[1])] > 0]
    log_fig = px.histogram(log2ipin_dataframe, x="Gene Name", y="Log2FC_{}_{}_IP_IN".format(group_list[0],
                                                                                            group_list[1]))
    log_fig.update_layout(xaxis=dict(showgrid=False, showticklabels=True, showline=True, linecolor='black'),
                          yaxis=dict(showgrid=False, showticklabels=True, showline=True, linecolor='black'),
                          plot_bgcolor='white')
    counter = counter + 1
    log_fig.write_html(filepath + os.sep + "{}_{}_{}_barplot_gene_name.html".format(counter, group_list[1],
                                                                                    group_list[0]))
    log_fig.write_image(filepath + os.sep + "{}_{}_{}_barplot_gene_name.pdf".format(counter, group_list[1],
                                                                                    group_list[0]))
    my_logger.info("Barplot is saved to {}, only the top genes are displayed here".format(
        filepath + os.sep + "{}_{}_{}_barplot_gene_name.html".format(counter, group_list[1], group_list[0])))

    counter = counter + 1
    tpm_percentile_qn_IN_and_IP_dataframe.to_csv(filepath+os.sep+'{}_ip_df_after_filtering_DE_genename'.format(counter)+'.csv', index=False)
    my_logger.info("Data frame is saved to {} after setting the threshold's for BC, SSMD and Log2FC".format(
        filepath+os.sep+'{}_ip_df_after_filtering_de_genename'.format(counter)+'.csv'))

    print("\n<---------------- Plots after the TRAP Analysis --------------->\n")
    my_logger.info("<---- Plots after the TRAP Analysis ---->")
    # tpm_percentile_qn_IN_and_IP_dataframe.sort_values('ssmd_ip_sni_sham', ascending=False, inplace=True)

    # [Reg Plot] scatter plot -- post TRAP processing
    plots_list = [tpm_percentile_qn_IN_and_IP_dataframe[sample2_qn_columns],
                  tpm_percentile_qn_IN_and_IP_dataframe[sample1_qn_columns],
                  tpm_percentile_qn_IN_and_IP_dataframe[sample2_qn_columns + sample1_qn_columns]]
    counter = plot_scatter(plots_list, plot_name_list, n, counter, filepath)

    # Scatter plots for SSMD and Log2FC based on threshold/filter/statistic
    temp_plot_df = tpm_percentile_qn_IN_and_IP_dataframe_rev_de_filter
    tpm_percentile_qn_IN_and_IP_dataframe_rev_de_filter.loc[
        (tpm_percentile_qn_IN_and_IP_dataframe_rev_de_filter['BC_{}_{}_IP'.format(group_list[0], group_list[1])].abs() >= bc_var) &
        (tpm_percentile_qn_IN_and_IP_dataframe_rev_de_filter['ssmd_ip_{}_{}'.format(group_list[1], group_list[0])].abs() >= ssmd_var) &
        (tpm_percentile_qn_IN_and_IP_dataframe_rev_de_filter['Log2fc_{}_{}_IP'.format(group_list[1], group_list[0])].abs() >= log2fc_var), 'Gene Expression'] = 'Differentially Expressed Genes'
    tpm_percentile_qn_IN_and_IP_dataframe_rev_de_filter['Gene Expression'].fillna('Genes not Differentially Expressed', inplace=True)
    plotly_fig = px.scatter(tpm_percentile_qn_IN_and_IP_dataframe_rev_de_filter,
                            x='Log2fc_{}_{}_IP'.format(group_list[1], group_list[0]),
                            y='ssmd_ip_{}_{}'.format(group_list[1], group_list[0]),
                            hover_data=['Gene Name'], color='Gene Expression', color_discrete_sequence=["black", "red"])
    plotly_fig.update_traces(marker_size=4)
    plotly_fig.update_layout(title='Scatter plot between SSMD and Log2FC',
                             xaxis=dict(showgrid=False, showticklabels=True, showline=True, linecolor='black'),
                             yaxis=dict(showgrid=False, showticklabels=True, showline=True, linecolor='black'),
                             plot_bgcolor='white')
    counter = counter + 1
    plotly_fig.write_html(filepath + os.sep + '{}_scatter_expression_IP_'.format(counter) + 'tpm.html')
    plotly_fig.write_image(filepath + os.sep + '{}_scatter_expression_IP_'.format(counter) + 'tpm.pdf')
    my_logger.info("Scatter plot for differentially expressed genes is plotted here {}".format(
        filepath + os.sep + '{}_scatter_expression_IP_'.format(counter) + 'tpm.html'))

    # Scatter plot of Log-IN and combined [input_threshold_df, temp_plot_df]
    qn_in_df_group1 = input_threshold_df.loc[:, sample1_qn_columns[:n]]
    qn_in_df_group2 = input_threshold_df.loc[:, sample1_qn_columns[n:]]
    input_threshold_df['Log2fc_{}_{}_IN'.format(group_list[1], group_list[0])] = traplib.log2_fold_change(qn_in_df_group1,
                                                                                           qn_in_df_group2, averaging='median')
    input_threshold_df['BC_{}_{}_IN'.format(group_list[0], group_list[1])] = pd.DataFrame(traplib.bhatt_distance(qn_in_df_group1,
                                                                                                  qn_in_df_group2))

    input_threshold_df.loc[
        (input_threshold_df['BC_{}_{}_IN'.format(group_list[0], group_list[1])].abs() >= bc_var) &
        (input_threshold_df['ssmd_in_{}_{}'.format(group_list[1], group_list[0])].abs() >= ssmd_var) &
        (input_threshold_df['Log2fc_{}_{}_IN'.format(group_list[1], group_list[0])].abs() >= log2fc_var), 'Gene Expression IN'] = 'Differentially Expressed Genes'
    input_threshold_df['Gene Expression IN'].fillna('Genes not Differentially Expressed', inplace=True)

    plotly_fig = px.scatter(input_threshold_df, x='Log2fc_{}_{}_IN'.format(group_list[1], group_list[0]),
                            y='ssmd_in_{}_{}'.format(group_list[1], group_list[0]), hover_data=['Gene Name'],
                            color='Gene Expression IN', color_discrete_sequence=["black", "red"])
    plotly_fig.update_traces(marker_size=4)
    plotly_fig.update_layout(title='Scatter plot between SSMD and Log2FC',
                             xaxis=dict(showgrid=False, showticklabels=True, showline=True, linecolor='black'),
                             yaxis=dict(showgrid=False, showticklabels=True, showline=True, linecolor='black'),
                             plot_bgcolor='white')
    counter = counter + 1
    plotly_fig.write_html(filepath + os.sep + '{}_scatter_expression_IN_'.format(counter) + 'tpm.html')
    plotly_fig.write_image(filepath + os.sep + '{}_scatter_expression_IN_'.format(counter) + 'tpm.pdf')
    my_logger.info("Scatter plot for differentially expressed genes is plotted here {}".format(
        filepath + os.sep + '{}_scatter_expression_IN_'.format(counter) + 'tpm.html'))

    # bar plot -- IP
    counter = plot_bar_chart(tpm_percentile_qn_IN_and_IP_dataframe, group_list, 'IP', counter, filepath)

    # bar plot -- IN
    input_threshold_df = input_threshold_df.loc[
                         (input_threshold_df['BC_{}_{}_IN'.format(group_list[0], group_list[1])].abs() >= bc_var) &
                         (input_threshold_df['ssmd_in_{}_{}'.format(group_list[1], group_list[0])].abs() >= ssmd_var) &
                         (input_threshold_df['Log2fc_{}_{}_IN'.format(group_list[1], group_list[0])].abs() >= np.log2(log2fc_var)), :]
    final_ip_df_bar_plot_data = input_threshold_df[['Gene Name', 'Log2fc_{}_{}_IN'.format(group_list[1], group_list[0])]]
    counter = counter + 1
    final_ip_df_bar_plot_data.to_csv(filepath + os.sep + "{}_percentile_threshold_IN_percentile_qn_tpm_Log2fcIN.csv".format(counter), index=False)

    counter = plot_bar_chart(final_ip_df_bar_plot_data, group_list, 'IN', counter, filepath)

    # generate heatmap post qn
    counter = counter + 1
    m = px.imshow(tpm_percentile_qn_IN_and_IP_dataframe[sample1_qn_columns+sample2_qn_columns].corr(),
                  color_continuous_scale="viridis", aspect=True)
    m.update_layout(title="Heatmap for {}".format('QN Corr Matrix'))
    m.write_html(filepath + os.sep + '{}_heatmap_qn'.format(counter) + '_tpm.html')
    m.write_image(filepath + os.sep + '{}_heatmap_qn'.format(counter) + '_tpm.pdf')
    my_logger.info("Heatmap for IP & IN is plotted here for the QN values {}".format(
        filepath + os.sep + '{}_heatmap_qn'.format(counter) + '_tpm.html'))

    # create z-score clustermap
    counter = counter + 1
    tpm_percentile_qn_IN_and_IP_dataframe = tpm_percentile_qn_IN_and_IP_dataframe[
        gene_name+sample1_qn_columns+sample2_qn_columns]
    tpm_percentile_qn_IN_and_IP_dataframe.set_index('Gene Name', inplace=True)
    de_in = sns.clustermap(tpm_percentile_qn_IN_and_IP_dataframe[sample1_qn_columns+sample2_qn_columns], cmap="YlGnBu",
                           metric="correlation", method="average", z_score=0, robust=False, figsize=(80, 80))
    de_in.savefig(filepath + os.sep + "{}_cmap_correlation_average_IP_{}_{}.pdf".format(counter, group_list[1], group_list[0]))
    my_logger.info("Cluster-map is plotted here {}".format(
        filepath + os.sep + '{}_heatmap_qn'.format(counter) + '_tpm.html'))

    my_logger.info("The analysis output is stored here {}".format(filepath))
    my_logger.info("<--------- The End --------->")
    print("The analysis output is stored here {}".format(filepath))
    print("<---------- The End ---------->")


# genes of interest
gene_id_list = ['Mrgprd', 'Calcb', 'Scn10a', 'Calca', 'Trpv1', 'Prph', 'Prdm12', 'Cx3cr1', 'Aif1', 'Aldh1l1', 'Cldn5', 'Mpz', 'Mbp', 'Gfap']
sample_name_list = ['IP_values', 'IN_values', 'IN_and_IP_values']

# #### EXAMPLES ##### group1-in and group2-ip
group1_in_list = ['77-10-IN-F259',	'77-14-IN-F26',	'77-3-IN-F081418',	'77-7-IN-F083018']
group1_ip_list = ['77-9-IP-F259',	'77-13-IP-F26',	'77-4-IP-F081418',	'77-8-IP-F083018']
group2_in_list = ['77-12-IN-M259',	'77-16-IN-M26',	'77-5-IN-M083018',	'77-1-IN-M081418']
group2_ip_list = ['77-11-IP-M259',	'77-15-IP-M26',	'77-6-IP-M083018',	'77-2-IP-M081418']

group_name_list = ['Exp-Group', 'Control-Group']  # [ group2, group1 ] or [ Experimental Group, Control Group ]

tpm_filepath = r"C:\Users\NXI220005\Downloads\new_trap_testing\TRAP-Analysis\trap_tpm_original.xlsx"

if __name__ == '__main__':

    """
    Note:
        Group 1 -- Experimental Group [Ex. SNI]
        Group 2 -- Control Group     [Ex. SHAM]
    """
    filepath = os.path.dirname(tpm_filepath) + os.sep + "Processed_files"
    my_logger = logging_functionality(filepath=filepath)

    # The code inside this main block will only execute if this script is run directly, not when imported as a module
    trap_analysis(group1_in_list, group1_ip_list, group2_in_list, group2_ip_list,
                  tpm_filepath,  # input file
                  gene_id_list, sample_name_list, group_name_list)
