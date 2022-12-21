import matplotlib.pyplot as plt

# uma dezena
alpha = 10
beta = 12
gamma = 14
# algumas dezenas
Ca = 20
Cb = 0
Cc = 25
# algumas unidades
Tmax = 1
T = 0
deltaT = 0.00025


def derivada_Ca(Ca, Cb, Cc, alpha):
    return -alpha * Ca * Cc + Cb


def derivada_Cb(Ca, Cb, Cc, beta):
    return beta * Ca * Cc - Cb


def derivada_Cc(Ca, Cb, Cc, gamma):
    return -gamma * Ca * Cc + Cb - 2 * Cc


def rungeKutta(T, deltaT, Ca, Cb, Cc, alpha, beta, gamma):
    k1 = [0, 0, 0]
    k2 = [0, 0, 0]
    k3 = [0, 0, 0]
    y = [Ca, Cb, Cc]
    vetor_Ca = [Ca]
    vetor_Cb = [Cb]
    vetor_Cc = [Cc]
    vetor_T = [T]

    while T < Tmax:

        k1[0] = derivada_Ca(Ca, Cb, Cc, alpha)
        k1[1] = derivada_Cb(Ca, Cb, Cc, beta)
        k1[2] = derivada_Cc(Ca, Cb, Cc, gamma)

        k2[0] = derivada_Ca(
            Ca + k1[0] * (deltaT / 2),
            Cb + k1[1] * (deltaT / 2),
            Cc + k1[2] * (deltaT / 2),
            alpha,
        )
        k2[1] = derivada_Cb(
            Ca + k1[0] * (deltaT / 2),
            Cb + k1[1] * (deltaT / 2),
            Cc + k1[2] * (deltaT / 2),
            beta,
        )
        k2[2] = derivada_Cc(
            Ca + k1[0] * (deltaT / 2),
            Cb + k1[1] * (deltaT / 2),
            Cc + k1[2] * (deltaT / 2),
            gamma,
        )

        k3[0] = derivada_Ca(
            Ca - (k1[0] * deltaT) + (2 * k2[0] * deltaT),
            Cb - (k1[1] * deltaT) + (2 * k2[1] * deltaT),
            Cc - (k1[2] * deltaT) + (2 * k2[2] * deltaT),
            alpha,
        )
        k3[1] = derivada_Cb(
            Ca - (k1[0] * deltaT) + (2 * k2[0] * deltaT),
            Cb - (k1[1] * deltaT) + (2 * k2[1] * deltaT),
            Cc - (k1[2] * deltaT) + (2 * k2[2] * deltaT),
            beta,
        )
        k3[2] = derivada_Cc(
            Ca - (k1[0] * deltaT) + (2 * k2[0] * deltaT),
            Cb - (k1[1] * deltaT) + (2 * k2[1] * deltaT),
            Cc - (k1[2] * deltaT) + (2 * k2[2] * deltaT),
            gamma,
        )

        y[0] = y[0] + (1 / 6) * (k1[0] + 4 * k2[0] + k3[0]) * deltaT
        y[1] = y[1] + (1 / 6) * (k1[1] + 4 * k2[1] + k3[1]) * deltaT
        y[2] = y[2] + (1 / 6) * (k1[2] + 4 * k2[2] + k3[2]) * deltaT

        Ca = y[0]
        Cb = y[1]
        Cc = y[2]

        vetor_Ca.append(Ca)
        vetor_Cb.append(Cb)
        vetor_Cc.append(Cc)
        T = T + deltaT
        vetor_T.append(T)

    return Ca, Cb, Cc, vetor_Ca, vetor_Cb, vetor_Cc, vetor_T


resultados_RK = rungeKutta(T, deltaT, Ca, Cb, Cc, alpha, beta, gamma)

print("\nDEPOIS DAS ITERAÇOES, AS APROXIMACOES FORAM:")
print(
    "\nCa = {}\nCb = {}\nCc = {}".format(
        resultados_RK[0], resultados_RK[1], resultados_RK[2]
    )
)


def grafico_RK(resultados_RK):

    plt.title("Gráfico dos valores de Ca,Cb,Cc com variação do tempo")
    plt.plot(resultados_RK[6], resultados_RK[3], color="blue", label="Ca")
    plt.plot(resultados_RK[6], resultados_RK[4], color="red", label="Cb")
    plt.plot(resultados_RK[6], resultados_RK[5], color="purple", label="Cc")
    plt.legend()
    plt.xlabel("Tempo(s)")
    plt.ylabel("Valor de C")
    plt.show()


grafico_RK(resultados_RK)


def refinamento():
    vetor_deltaT = [deltaT]

    for i in range(0, 4):
        vetor_deltaT.append(vetor_deltaT[i] * 2)

    fig, (grafico_Ca, grafico_Cb, grafico_Cc) = plt.subplots(3, 1)

    for i in range(0, 5):
        resultados_RK = rungeKutta(T, vetor_deltaT[i], Ca, Cb, Cc, alpha, beta, gamma)
        cor = ["blue", "pink", "purple", "red", "black"]

        grafico_Ca.plot(
            resultados_RK[6],
            resultados_RK[3],
            color=cor[i],
            label="%.4f" % (float(vetor_deltaT[i])),
        )
        grafico_Ca.legend(fontsize=7)
        grafico_Ca.set_xlabel("tempo(s)", fontsize=9)
        grafico_Ca.set_ylabel("Valor de Ca", fontsize=9)

        grafico_Cb.plot(
            resultados_RK[6],
            resultados_RK[4],
            color=cor[i],
            label="%.4f" % (float(vetor_deltaT[i])),
        )
        grafico_Cb.legend(fontsize=7)
        grafico_Cb.set_xlabel("tempo(s)", fontsize=9)
        grafico_Cb.set_ylabel("Valor de Cb", fontsize=9)

        grafico_Cc.plot(
            resultados_RK[6],
            resultados_RK[5],
            color=cor[i],
            label="%.4f" % (float(vetor_deltaT[i])),
        )
        grafico_Cc.legend(fontsize=7)
        grafico_Cc.set_xlabel("tempo(s)", fontsize=9)
        grafico_Cc.set_ylabel("Valor de Cc", fontsize=9)

    plt.suptitle("Gráfico de refinamento para diferentes valores de deltaT")
    plt.show()


refinamento()


def sensibilidade():

    sensibilidade_alpha = [alpha]
    sensibilidade_beta = [beta]
    sensibilidade_gamma = [gamma]

    i = 0
    for i in range(0, 4):
        sensibilidade_alpha.append(sensibilidade_alpha[i] / 2)
        sensibilidade_beta.append(sensibilidade_beta[i] / 2)
        sensibilidade_gamma.append(sensibilidade_gamma[i] / 2)

    fig, (grafico_alpha, grafico_beta, grafico_gamma) = plt.subplots(3, 1)

    for i in range(0, 5):
        resultados_RK = rungeKutta(
            T,
            deltaT,
            Ca,
            Cb,
            Cc,
            sensibilidade_alpha[i],
            sensibilidade_beta[i],
            sensibilidade_gamma[i],
        )
        cor = ["blue", "pink", "purple", "red", "black"]

        grafico_alpha.plot(
            resultados_RK[6],
            resultados_RK[3],
            color=cor[i],
            label="%.2f" % (sensibilidade_alpha[i]),
        )
        grafico_alpha.legend(title="alpha", fontsize=7)
        grafico_alpha.set_xlabel("tempo(s)", fontsize=9)
        grafico_alpha.set_ylabel("Valor de Ca", fontsize=9)

        grafico_beta.plot(
            resultados_RK[6],
            resultados_RK[4],
            color=cor[i],
            label="%.2f" % (sensibilidade_beta[i]),
        )
        grafico_beta.legend(title="beta", fontsize=7)
        grafico_beta.set_xlabel("tempo(s)", fontsize=9)
        grafico_beta.set_ylabel("Valor de Cb", fontsize=9)

        grafico_gamma.plot(
            resultados_RK[6],
            resultados_RK[5],
            color=cor[i],
            label="%.2f" % (sensibilidade_gamma[i]),
        )
        grafico_gamma.legend(title="gamma", fontsize=7)
        grafico_gamma.set_xlabel("tempo(s)", fontsize=9)
        grafico_gamma.set_ylabel("Valor de Cc", fontsize=9)

    plt.suptitle(
        "Gráfico de sensibilidade para diferentes valores de alpha, beta e gamma"
    )
    plt.show()


sensibilidade()
