import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

sns.set_style(style="whitegrid")

data_file = snakemake.input.full_sno_host

df = pd.read_csv(data_file, sep="\t")

print(df.columns)

df_gb = df[['merged_name', 'interaction_type']].groupby('interaction_type').count().reset_index()
order = {
    'intra': 0,
    'next_exon': 1,
    'previous_exon': 2,
    'distant': 3
}
df_gb['order'] = df_gb.interaction_type.map(order)
df_gb.sort_values('order', inplace=True)

data = dict(zip(df_gb.interaction_type, df_gb.merged_name))


print(df_gb)


fig, ax= plt.subplots(figsize=(8,8))
ax.axis('equal')

new_names = {
    'intra': 'Same intron',
    'next_exon': 'Next exon',
    'previous_exon': 'Previous exon',
    'distant': 'Distant'
}
labels = [new_names[x] for x in data.keys()]
values = data.values()
colors = ['#b2182b', '#ef8a62', '#fddbc7', '#e0e0e0']

# Wedge properties
wp = { 'linewidth' : 2, 'edgecolor' : "white" }

def autopct(val):
    num = (val / 100) * sum(values)
    # return f'{val:.1f}% ({num:.0f})'
    return f'{val:.1f}%'

wedges, texts, autotexts = ax.pie(values, labels=labels,
       autopct=autopct, colors=colors,
       shadow=False, startangle=0,
       textprops={'color':'black', 'fontsize': 20},
       wedgeprops=wp)

ax.set_title('Proportion of location of snoRNA targets', fontsize=25)

ax.legend(wedges, labels,
      title ="snoRNA target location",
      loc="upper left",
      bbox_to_anchor=(1, 0, 1, 1),
      fontsize=16)

fig.subplots_adjust(right=0.7)

# plt.show()
plt.savefig(snakemake.output.svg, format='svg')
