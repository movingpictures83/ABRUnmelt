import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("whitegrid")
plt.ylim([0,20])

import PyPluMA
import PyIO

class ABRUnmeltPlugin:
  def input(self, inputfile):
      self.parameters = PyIO.readParameters(inputfile)

  def run(self):
      pass

  def output(self, outputfile):
      amr_counts_file = PyPluMA.prefix()+"/"+self.parameters["amrcounts"]#'amr_counts_metagen/amr_counts_merged.csv'
      read_counts_file = PyPluMA.prefix()+"/"+self.parameters["readcounts"]#'read_counts.csv'
      metadata_file = PyPluMA.prefix()+"/"+self.parameters["metadata"]#'read_counts.csv'
      read_counts = pd.read_csv(read_counts_file)
      read_counts["Samples#"] = read_counts["sample"]
      read_counts = read_counts[['Samples#', 'n_reads']]

      amr_counts = pd.read_csv(amr_counts_file)
      amr_counts = amr_counts.merge(read_counts, on="Samples#")

      amr_counts['RPKM'] = (amr_counts['read_count'] * 1000 * 1000000)/(amr_counts['n_reads'] * amr_counts['gene_length'])

      df_unmelted = amr_counts[['Samples#', 'RPKM', 'gene_name', "Group"]]

      # Unmelt with gene granularity
      #df_unmelted = df_unmelted.set_index(['sampleID', 'diagnosis', 't(timepoint)'])['RPKM'].unstack().reset_index()
      print(amr_counts)
      df_unmelted = df_unmelted.pivot_table(index=['Samples#', 'Group'], columns='gene_name', values='RPKM')
      df_unmelted= df_unmelted.fillna(0)

      df_unmelted = df_unmelted.reset_index()

      # Add samples with zero profile
      metadata_df = pd.read_csv('samples_metagen_metadata.csv')
      metadata_df = metadata_df[['Samples#', 'Group']]
      df_unmelted = metadata_df.merge(df_unmelted, how='left')
      df_unmelted = df_unmelted.fillna(0)

      df_unmelted.to_csv(outputfile+"_unmelted.csv")

      print(metadata_df)


      # melt the dataframe
      melted_df = pd.melt(df_unmelted, id_vars=['Samples#', 'Group'])

      melted_df['RPKM'] = melted_df['value']
      melted_df['gene_name'] = melted_df['variable']
      melted_df = melted_df.drop('variable', axis=1)
      melted_df = melted_df.drop('value', axis=1)
      melted_df.to_csv(outputfile+"_counts_metagen_with_zeros.csv", index=False)

      #sns.boxplot(x='gene_name', y='RPKM', hue='Group', data=melted_df)
      #plt.ylim(0,15)
      #plt.show()

