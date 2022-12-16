# uma dezena
beta = 1
alpha = 1
gama = 1
# algumas dezenas
Ca = 1
Cb = 0
Cc = 1
# algumas unidades
Tmax = 0.5
T = 0
deltaT = 0.5

k1 = [0, 0, 0]
k2 = [0, 0, 0]
k3 = [0, 0, 0]
y = [Ca, Cb, Cc]

derivada_Ca = -alpha*Ca*Cc + Cb #1funçao
derivada_Cb = beta*Ca*Cc - Cb #2funçao
derivada_Cc = -gama*Ca*Cc + Cb - 2*Cc #3funçao

def derivada_Ca(Ca, Cb, Cc, alpha): 
    return (-alpha*Ca*Cc + Cb) 

def derivada_Cb(Ca, Cb, Cc, beta): 
    return (beta*Ca*Cc - Cb) 

def derivada_Cc(Ca, Cb, Cc, gama): 
    return (-gama*Ca*Cc + Cb - 2*Cc) 

def rungeKutta(T, deltaT, Ca, Cb, Cc, alpha, beta, gama): 

    while(T<=Tmax):

        print("\nITERAÇÃO DO TEMPO {}".format(T))

        k1[0] = derivada_Ca(Ca, Cb, Cc, alpha)
        k1[1] = derivada_Cb(Ca, Cb, Cc, beta)
        k1[2] = derivada_Cc(Ca, Cb, Cc, gama)
        print("k1_Ca = {}\nk1_Cb = {}\nk1_Cc = {}".format(k1[0], k1[1], k1[2]))

        k2[0] = derivada_Ca(Ca+k1[0]*(deltaT/2), Cb+k1[1]*(deltaT/2), Cc+k1[2]*(deltaT/2), alpha)
        k2[1] = derivada_Cb(Ca+k1[0]*(deltaT/2), Cb+k1[1]*(deltaT/2), Cc+k1[2]*(deltaT/2), beta)
        k2[2] = derivada_Cc(Ca+k1[0]*(deltaT/2), Cb+k1[1]*(deltaT/2), Cc+k1[2]*(deltaT/2), gama)
        print("\nk2_Ca = {}\nk2_Cb = {}\nk2_Cc = {}".format(k2[0], k2[1], k2[2]))


        k3[0] = derivada_Ca(Ca-(k1[0]*deltaT)+(2*k2[0]*deltaT), Cb-(k1[1]*deltaT)+(2*k2[1]*deltaT), Cc-(k1[2]*deltaT)+(2*k2[2]*deltaT), alpha)
        k3[1] = derivada_Cb(Ca-(k1[0]*deltaT)+(2*k2[0]*deltaT), Cb-(k1[1]*deltaT)+(2*k2[1]*deltaT), Cc-(k1[2]*deltaT)+(2*k2[2]*deltaT), beta)
        k3[2] = derivada_Cc(Ca-(k1[0]*deltaT)+(2*k2[0]*deltaT), Cb-(k1[1]*deltaT)+(2*k2[1]*deltaT), Cc-(k1[2]*deltaT)+(2*k2[2]*deltaT), gama)
        print("\nk3_Ca = {}\nk3_Cb = {}\nk3_Cc = {}".format(k3[0], k3[1], k3[2]))


        y[0] = y[0] + (1/6)*(k1[0]+4*k2[0]+k3[0])*deltaT
        y[1] = y[1] + (1/6)*(k1[1]+4*k2[1]+k3[1])*deltaT
        y[2] = y[2] + (1/6)*(k1[2]+4*k2[2]+k3[2])*deltaT
        print("\nCa = {}\nCb = {}\nCc = {}".format(y[0], y[1], y[2]))

        Ca = y[0]
        Cb = y[1]
        Cc = y[2]

        T = T+deltaT


    return y

rungeKutta(T, deltaT, Ca, Cb, Cc, alpha, beta, gama)
