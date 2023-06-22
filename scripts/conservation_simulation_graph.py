import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scipy.stats as stats

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 15})
plt.rcParams['font.sans-serif'] = ['Arial']

def fisher_exact(net, simulated):

    odds_ratios = []
    pvalues = []
    for sim in simulated:
        oddsratio, pvalue = stats.fisher_exact([[*net[0]],
                                                [*sim]])
        print(oddsratio, pvalue)
        odds_ratios.append(oddsratio)
        pvalues.append(pvalue)

    print('---------------------------')
    print(f'ods ratio: {np.mean(odds_ratios)}, pValue: {np.mean(pvalues)}')


def graph():

    def percent(list_):
        return np.array([
            [(x/(x+y))*100, (y/(x+y))*100]
            for x, y in list_
        ])


    # network_val_ = [(8, 68)] # without duplicates
    network_val_ = [(11, 91)] # with all
    network_val = percent(network_val_)
    # Simulated with 0.5
    simulated_vals_ = [
        (19, 433),
        (17, 435),
        (14, 438),
        (15, 437),
        (14, 438),
        (15, 437),
        (13, 439),
        (13, 439),
        (15, 437),
        (18, 434)
    ]
    # --------------------------------------------
    # Simulated with 0.2
    # network_val_ = [(16, 86)]
    # simulated_vals_ = [
    #     (57, 395),
    #     (47, 405),
    #     (59, 393),
    #     (62, 390),
    #     (51, 401),
    #     (57, 395),
    #     (48, 404),
    #     (57, 395),
    #     (56, 396),
    #     (63, 389)
    # ]
    # --------------------------------------------
    simulated_vals = percent(simulated_vals_)

    lows = [network_val[0, 1], np.mean(simulated_vals[:,1])]
    highs = [network_val[0, 0], np.mean(simulated_vals[:,0])]

    THRESH = 0.2
    N = 2
    ind = np.arange(N)    # the x locations for the groups
    width = 0.75       # the width of the bars: can also be len(x) sequence

    fig1, ax = plt.subplots(figsize=(8, 8))
    # ax.bar(ind, lows, width,label=f'mean cons >= {THRESH}', color='#377eb8')
    # ax.bar(ind, highs, width, bottom=lows, label=f'mean cons >= {THRESH}', color='#e41a1c')
    ax.bar(ind, highs, width, label=f'mean cons >= {THRESH}', color='#80b1d3',
            yerr=[0,.5], error_kw={'capsize' : 7})


    plt.ylabel('Proportion with high mean conservation (%)')
    plt.title('Mean conservation of snoRNA target regions')
    plt.xticks(ind, ('snoRNA interacting\nregion in host intron', 'Random'))#, rotation='vertical')
    # plt.subplots_adjust(bottom=0.2, right=0.8, left=0.08)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.savefig(snakemake.output.svg, format='svg')
    # plt.show()

    fisher_exact(network_val_, simulated_vals_)


def main():

    graph()


if __name__ == '__main__':
    main()
