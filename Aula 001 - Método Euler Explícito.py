

"""
OBS: Não se assuste com a ordem das aulas, por exemplo, Aula 7 e depois Aula 13. O que acontece é que temos aulas de
     dúvidas e aulas de exercícios, que não colocarei aqui como aula. Terá uma pasta chamada Lista 1,2,3... que você
     podera encontrar os exercícios para treinar.
Aula 1
A equação abaixo é um exemplo de equação diferêncial:
    d²x/d²t+q(t)/dt dx = r(x)
- Onde t é a variável independente e x(t) é a var. dependente.
- Poderia ser uma força dissipativa entre outros vários tipos.
- São sempre em relação à uma variável independente (neste caso o tempo).
Podemos pensar na condição de contorno da equação como se fosse uma barra com calor
em ambas pontas, porém, não trabalharemos com este tipo de equação, mas sim com equações
com valores iniciais.
Solução Numérica.
    Método Euler (não deve ser utilizado).
        - Transformar um dx > delta_x = x(n+1) + 1 - xn
        - Transformar um dt > delta_t = t(n+1) + 1 - tn
    Obtemos:
        - detlta_x/delta_x = f(t; x)
    Usando a definição:
        - x(n+1) = xn + delta_t f(t; x)
# Método de Euler Explícito: (não deve ser utilizado no dia a dia)
Prob. inicial.
dx/dt = x
x(0) = x0
delta_x/delta_t = x = f(t, x) = f(x)
x(n+1) = xn + xn * delta_t
Onde x*delta_t é o Método de Euler Explícito.
Ficamos com:
    x(n+1) = (1 + delta_t)*xn
"""

# Exemplo Teste:

"""
Fórmula da dilatação térmica:
 
  delta_L = L0 * alfa * delta_T
ou
  L = L0 + L0 * alfa * delta_T
  
Onde dL/dT = alfa * L
Exercício: No que resulta mantermos T0 e Tf fixos e aumentarmos o número de
passos entre estes dois valores?
L0 = 1, alfa = 0.01, T0 = 10 , Tf = 50
"""
import numpy as np
import matplotlib.pyplot as plt


def dilatacao_termica(n, L0=1, alfa=0.01, T0=10, Tf=50):
    """
    Obejtivo é utilizar 'n' para o número de passos e ver a diferença no
    resultado final, junto de seu erro absoluto e relativo.
    """
    comprimento = [L0]

    delta_T = np.linspace(10, 50, num=n + 1)

    for i in range(0, n):
        L = L0 + L0 * alfa * (delta_T[i + 1] - delta_T[i])
        L0 = L
        comprimento.append(L)

    return comprimento


print(dilatacao_termica(1))
print(np.linspace(0, 60, 2))

plt.plot(np.linspace(10, 50, num=2), dilatacao_termica(1), label='1 Passo', color='crimson')
plt.plot(np.linspace(10, 50, num=11), dilatacao_termica(10), label='10 Passo', color='orange')
plt.plot(np.linspace(10, 50, num=101), dilatacao_termica(100), label='100 Passo', color='gold')
plt.xlabel('Temperatura')
plt.ylabel('Comprimento')
plt.title('Dilatação Térmica')
plt.legend()
plt.show()