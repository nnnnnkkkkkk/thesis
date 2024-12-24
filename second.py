import first as f
from numpy import roots, conj, log
from collections import Counter
import time
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool
import csv

# plnms = f.polynomialsOfDegreeN(3)
# df = pd.read_csv('/Users/nikolmetry/PycharmProjects/cusra/rootsP0degree4.csv')
# pls, rs, ab = [], [[], [], []], []
# # # shuffle(plnms)
# for g in plnms:
#     g = f.polynomial(g)
#     if g.isMonic() and g not in refplnms:
#         print(g)
#         # a_n is a coefficient in front of u^n in P_i
#         a0 = f.G(g, 1)
#         a3 = a0 * (- f.q ** 4 + f.q ** 3) + f.G(g, 4)
#         # sanity check, see p.22 and p.24 for the reference
#         c1 = conj(f.shiftedGaussSum(f.polynomial([1]), g)) * f.q ** (g.degree() * (-2 / 3)) * \
#              f.rho((1 - 2 * g.degree()) % 3)
#         c2 = a3 * f.q ** (-4) + a0
#         if abs(c1 - c2) > 1e-10:
#             print('!!!', g)
#             print(c1, c2)
#             if f.gcd(g, g.derivative()) == 1:
#                 break
#         else:
#             if abs(c1) > 1e-10:
#                 pls.append(g)
#                 print(c1, c2)
#                 r = sorted(roots([a3, 0, 0, a0]))
#                 for i in range(len(rs)):
#                     rs[i].append(r[i])
#                 ab.append(abs(r[0]))
#         # if len(pls) == 100:
#         #     break
#         # else:
#         #     if f.gcd(g, g.derivative()) == 1:
#         #         pls.append(g)
#         #         r = sorted(roots([a3, 0, 0, 1]))
#         #         for i in range(len(rs)):
#         #             rs[i].append(r[i])
#         #         ab.append(abs(r[0]))
#
# data = {"polynomial": pls, "1st root": rs[0], "2nd root": rs[1], "3rd root": rs[2], "root size": ab}
# frame = pd.DataFrame(data)
# print(frame)
# frame.to_csv("/Users/nikolmetry/rootsP1degree3.csv")
#
#
#
# df = df.sort_values(by=['root size'])
# fig, ax = plt.subplots()
# x, y, z, root_size = [], [], [], []
# curr = 0
# for i in df.index:
#     if abs(curr - df['root size'][i]) > 1e-11:
#         # print(curr, df['root size'][i])
#         root_size.append((df['root size'][i]))
#         x.append([])
#         y.append([])
#         curr = df['root size'][i]
#     x[int(len(root_size)) - 1].extend(
#         [complex(df['1st root'][i]).real, complex(df['2nd root'][i]).real, complex(df['3rd root'][i]).real])
#     y[int(len(root_size)) - 1].extend(
#         [complex(df['1st root'][i]).imag, complex(df['2nd root'][i]).imag, complex(df['3rd root'][i]).imag])
#     z.append(str(log(df["root size"][i])/log(7))[:12])
#
# labels, values = zip(*sorted(Counter(z).items()))
# indexes = arange(len(labels))
# width = 1
# labels = list(labels)
# for i in range(len(labels)):
#     labels[i] = labels[i][:8]


# for i in range(len(root_size)):
#     plt.scatter(x[i], y[i])
#
# plt.scatter(0,0, color='black')
# plt.annotate('(0, 0)', (0.014, 0.01))
# c1 = plt.Circle((0, 0), 1 / f.q, fill=False)
# c2 = plt.Circle((0, 0), 7 ** (-4 / 3), fill=False)
# c3 = plt.Circle((0, 0), 7 ** (-3 / 2), fill=False)
# c4 = plt.Circle((0, 0), 7 ** (-2), fill=False)
# ax.add_patch(c1)
# ax.add_patch(c2)
# ax.add_patch(c3)
# ax.add_patch(c4)
# ax.axvline(x=0, color='black')
# ax.axhline(y=0, color='black')
# ax.set_xlabel('x')
# ax.set_ylabel('y')


#
# ax.bar(indexes, values, width)
# ax.set_xticks(indexes, labels)
# ax.set_xlabel('Exponent of a root size')
# ax.set_ylabel('Number of roots')
# print((Counter(z)))
# plt.show()

if __name__ == "__main__":
    procs = 6
    # df = pd.read_csv('/Users/nikolmetry/PycharmProjects/cusra/rootsP0degree4.csv')
    # sqfreem = df['polynomial'][8:10]
    # i = 1
    # for g in sqfreem:
    #     g = f.polynomial(g)
    #     fs = f.polynomialsOfDegreeN(5)[:7 ** 5]
    #     shs = [g for _ in range(len(fs))]
    #     args = list(zip(shs, fs))
    #     start = time.time()
    #     pool = Pool(procs).starmap(f.monGaussSum, args)
    #     s = sum(pool)
    #     end = time.time()
    #     print(g, s, end - start)
    #     a0 = f.G(g, 2)
    #     a3 = a0 * (- f.q ** 4 + f.q ** 3) + s
    #     c2 = a3 * f.q ** (-4) + a0
    #     c1 = conj(f.shiftedGaussSum(f.polynomial([1]), g)) * f.q ** (
    #             g.degree() * (-2 / 3)) * \
    #          f.rho((2 - 2 * g.degree()) % 3) * 7 ** (8 / 3)
    #     if abs(c1 - c2) > 1e-10:
    #         print('!!!', g)
    #         print(c1, c2)
    #         if f.gcd(g, g.derivative()) == 1:
    #             break
    #     r = sorted(roots([a3, 0, 0, a0]))
    #     fields = [8 + i, g, r[0], r[1], r[2], abs(r[2])]
    #     with open('/Users/nikolmetry/PycharmProjects/cusra/rootsP2degree4.csv', 'a') as file:
    #         writer = csv.writer(file)
    #         writer.writerow(fields)
    #     i += 1

    g = f.polynomial([1, 0, 0, 2, 3])
    fs = f.polynomialsOfDegreeN(5)[:7 ** 5]
    shs = [g for _ in range(len(fs))]
    args = list(zip(shs, fs))
    start = time.time()
    pool = Pool(procs).starmap(f.monGaussSumNew, args)
    s = sum(pool)
    end = time.time()
    print(g, s, end - start)
    a0 = f.G(g, 2)
    a3 = a0 * (- f.q ** 4 + f.q ** 3) + s
    c2 = a3 * f.q ** (-4) + a0
    c1 = conj(f.shiftedGaussSum(f.polynomial([1]), g)) * f.q ** (
            g.degree() * (-2 / 3)) * f.q ** (4 * 2 / 3) * f.rho((2 - 2 * g.degree()) % 3)
    if abs(c1 - c2) > 1e-10:
        print('!!!', g)
        print(c1, c2)
        print(f.gcd(g, g.derivative()))
    else:
        r = sorted(roots([a3, 0, 0, a0]))
        fields = [0, g, r[0], r[1], r[2], abs(r[2])]
        with open('/Users/nikolmetry/PycharmProjects/cusra/rootsP2degree4.csv', 'a') as file:
            writer = csv.writer(file)
            writer.writerow(fields)
