# ---------------------------------------------------------------------------------------------
# Dados de entrada:
# ---------------------------------------------------------------------------------------------
av = 0  # usar pressão do reservatório na fronteira externa
beta = 0  # Coeficiente de Forchheimer, m^-1. Não relevante para o fluido 1.
f_w_sc = 0.0  # Corte de água, não relevante para fluidos 1 e 2
fluid = 3  # Tipo de fluido: 1 = óleo monofásico, 2 = gás monofásico, >=3 = fluxo multifásico
h = 20.0  # Altura do reservatório, m
k = 1.11e-13  # Permeabilidade efetiva para fluidos 1 e 2, ou absoluta para fluido 3, m^2
n_pt = 100  # Número de pontos no gráfico
oil = 1  # Modelo de óleo: 1 = black oil (Standing), 2 = black oil (Glaso), 3 = óleo volátil
p_R = 30.0e6  # Pressão do reservatório, Pa
q_o_sc_max = -5e-2  # Vazão máxima de óleo, m^3/s
r_e = 500  # Raio externo, m
r_w = 0.2  # Raio do poço, m
R_go = 200  # Razão gás-óleo, m^3/m^3. Não relevante para fluido 2
rho_g_sc = 0.95  # Densidade do gás nas condições padrão, kg/m^3
rho_o_sc = 850  # Densidade do óleo nas condições padrão, kg/m^3
rho_w_sc = 1000  # Densidade da água nas condições padrão, kg/m^3
semi = 0  # 0 = estado estacionário, 1 = estado semi-estacionário
simp = 1  # 0 = solução numérica, 1 = solução simplificada (semi-analítica)
S = 0  # Fator de dano (skin)
T_R = 60  # Temperatura do reservatório, °C

# Dados de permeabilidade relativa (não relevantes para fluidos 1 e 2)
k_rg_0 = 0.7  # Permeabilidade relativa do gás
k_ro_0 = 0.9  # Permeabilidade relativa do óleo
k_rw_0 = 0.5  # Permeabilidade relativa da água
n_g = 3  # Expoente de Corey para gás
n_og = 3  # Expoente de Corey para óleo em fluxo gás-óleo
n_ow = 3  # Expoente de Corey para óleo em fluxo óleo-água
n_w = 3  # Expoente de Corey para água
S_gc = 0.00  # Saturação crítica de gás
S_or = 0.10  # Saturação residual de óleo
S_wi = 0.15  # Saturação inicial de água

# Verificação de entrada
def verificar_entrada():
    global simp, semi, av
    if simp == 1 and beta != 0:
        print("Aviso: Fluxo de Forchheimer disponível apenas numericamente.")
        simp = 0
    if simp == 0 and semi == 1:
        print("Aviso: Solução de estado semi-estacionário disponível apenas analiticamente.")
        semi = 0
    if simp == 0 and av == 1:
        print("Aviso: Solução com pressão média do reservatório disponível apenas analiticamente.")
        av = 0

verificar_entrada()

# Criar vetores de dados
rel = [k_rg_0, k_ro_0, k_rw_0, n_g, n_og, n_ow, n_w, S_gc, S_or, S_wi]
rho_sc = [rho_g_sc, rho_o_sc, rho_w_sc]

# Função para calcular a curva IPR (exemplo simplificado)
def calcular_ipr():
    delta_q_o_sc = q_o_sc_max / n_pt  # Incremento de taxa de fluxo de óleo, m^3/s
    p_wf = p_R
    resultados = []
    
    for i in range(1, n_pt + 1):
        q_o_sc = i * delta_q_o_sc  # Taxa de óleo, m^3/s
        q_g_sc = R_go * q_o_sc  # Taxa de gás, m^3/s
        q_w_sc = (f_w_sc / (1 - f_w_sc)) * q_o_sc if f_w_sc > 0 else 0  # Taxa de água, m^3/s
        
        # Aqui deveria ser chamada a função res(), que precisa ser implementada
        # p_wf = res(beta, fluid, h, k, oil, p_R, [q_g_sc, q_o_sc, q_w_sc], r_e, r_w, rel, rho_sc, T_R)
        p_wf = p_R - i * 1e5  # Exemplo de cálculo fictício (substituir por função real)
        
        if p_wf < 1e5:  # Evitar pressões abaixo da atmosférica
            break
        resultados.append([-q_o_sc, p_wf])
    
    return np.array(resultados)

# Calcular e plotar a curva IPR
resultados = calcular_ipr()
if len(resultados) > 0:
    plt.plot(resultados[:, 0] * 1e3, resultados[:, 1] / 1e6, linewidth=1)
    plt.xlabel('Vazão de Óleo, -q_{o,sc} (10^{-3} m³/s)')
    plt.ylabel('Pressão de Fundo, p_{wf} (MPa)')
    plt.grid(True)
    plt.show()