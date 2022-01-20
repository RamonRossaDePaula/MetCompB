"""
Método de Euler Implícito.
- No método anterior encontramos o valor de x(n+1) utilizando o valor encontrado anteriormente ( x(n) ),
neste método agora encontraremos o x(n+1) a partir de f(t(n+1), x(n+1)), valor este que não foi encontrado,
por isso o método é IMPLÍCITO.

x(n+1) = xn + f( t(n+1), x(n+1) ) * h

• Ponto Favorável do Método:
    - A partir do mesmo valor de passo, delta_t = h, ele dará um valor da solução melhor do que Euler Explícito.

• Ponto Desfavorável do Método:
    - Por ser Implícito, não temos o termo x(n+1) isolado do lado esquerdo.
    - Dependendo da função f(x) (função seno, exponencial, tangente...) funções *i*transcendentais, pode ser que
      não tenha solução analítica. Esta equação pode ser calculada analiticamente ou, em casos mais complicados,
      numericamente através de um algoritmo de Newton-Raphson, por exemplo.
 *i*  Uma função transcendente (em inglês: transcendental) é uma função a qual não satisfaz uma equação polinomial
cujos coeficientes são eles próprios polinomiais. Em outras palavras, uma função é dita transcendente quando ela
não pode ser expressa por uma combinação finita de operações algébricas.
"""


"""
Exemplo:
Calcular o decaimento radiotivo dNu/dT = -Nu*t.
• Pseudocódigo 1;
    - Declarar as variáveis; Ni, N0, delta_t
    - Mudar variáveis; N = Ni, t = t0 = 0, delta_t = h
    - Fazer o cálculo.
    - Guardar os resultados.
• Pseudocódigo 2; Inicialização.
    - Dar valor às variáveis; N0 = 1000, tau = 0.1, h = 0.1
    - Colocar o tempo inicial t(0).
    - Colocar número de passos.
    
• Pseudocódigo 3; Calculando.
    - Para cado passo i (começando de i = 0), calcule Nu e t no passo i+1;
        - Nu(ti+1) = Nu(ti) - (Nu(ti)/tau)dt [Utilizar método de Euler]
        - ti+1 = ti + delta_t
        - Repita para n-1 passos
"""

"""
Exercício: Decaimento do Rubídio(Rb**82)
OBS: Os passos abaixo estão listando o Método Explícito.
dN/dt = - alfa * N
1 - Use;
        
        delta_N = - alfa * delta_t * N(t)
        delta_N = N(t + delta_t) - N(t) = - alfa * delta_t * N(t)
        N(t + delta_N) = (1 - alfa * delta_t) * N(t)
        
com delta_t = 1s, alfa = 9,24 * 10 ** -3 segundos...
...mostre que uma amostra de Rubídio terá sua radioatividade reduzida à metade depois de decorridos 75 s.
OBS 2: Refazer porém utilizando o método Implícito.
"""

import numpy as np
import matplotlib.pyplot as plt


def decaimento(n, t, N0 = 1000, alfa = 9.24e-3):

    metade = [N0]

    delta_t = np.linspace(0, t, num =n + 1)

    for i in range(0,n):
        N = N0 - alfa * N0 * (delta_t[i + 1] - delta_t[i])
        N0 = N
        metade.append(N)

    return metade

#decaimento(100,75)

"""
2 - Construa um gráfico com delta_t = 1 e N(0) = 1000 de N(t) em função do tempo.
Assinale no gráfico os valores obtidos a cada 10s.
"""

y = decaimento(74, 75)
x = np.linspace(0, 75, num =75)

pontos_x = []
pontos_y = []
for i in range(0,70,10):
    pontos_x.append(i)
    pontos_y.append(y[i])

plt.plot(x, y, color ='crimson')
plt.plot(pontos_x, pontos_y, color = 'black', marker = '.', linewidth=0)
plt.title('Decaimento do Rubídio')
plt.xlabel('t')
plt.ylabel('N(t)')
plt.grid(True)
plt.show()

"""
3 - Repita o gráfico com escala logarítmica na vertical.
"""

plt.plot(x, y, color ='crimson')
plt.plot(pontos_x, pontos_y, color = 'black', marker = '.', linewidth=0)
plt.title('Decaimento do Rubídio')
plt.xlabel('t')
plt.ylabel('N(t)')
plt.yscale('log')
plt.show()

"""
4 - Meça o coeficiente ângular da reta no caso anterior e compare com "alfa" do Rb.
"""