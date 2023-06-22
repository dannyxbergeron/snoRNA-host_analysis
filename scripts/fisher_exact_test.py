import numpy as np
import scipy.stats as stats


distances = [100, 200, 300, 500, 1000, 10000]
pos_net = [11, 15, 17, 20, 21, 29]
pos_others = [12, 22, 28, 40, 55, 115]

total_net = 60
total_other = 279

for dist, net, other in zip(distances, pos_net, pos_others):

    net_neg = total_net - net
    other_neg = total_other - other

    obs = np.array([[net, net_neg], [other, other_neg]])
    odds, p_val = stats.fisher_exact(obs)

    print(f'Distance from the splicing event: {dist}')
    print(f'Odds: {odds:.2f}, pValue: {p_val:.5f}')
    print('---------------------------------------------\n')
