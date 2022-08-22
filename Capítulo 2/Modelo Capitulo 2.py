from random import SystemRandom 
import pickle as pck 
import math
import os
random=SystemRandom()
######
#P - Q - N - H or Q - H
#P- parasitoide especialista
#Q- parasitoide generalista
#N- host atacado pelo P e pelo Q
#H- host atacado pelo Q
######
#a1 - eficiencia de busca do P pro N
#a2 - eficiencia de busca do Q pro N
#a3 - eficiencia de busca do Q pro H
#lambda1 - taxa de crescimento populacionao de N
#lambda2 - taxa de crescimento populacionao de H
#H0=é o numero de hosts que precisa ter para impedir das femeas de parasitoide se dispersarem ([PQ]) #to-do: esse valor precisa ter mais 1 pq pode ser diferente para o Q predando outro host?
#h0=é a tolerancia a outros hosts ([NH])
#f0=é a tolerancia a outros parasitoides femeas.([PQ]) #to-do: esse valor precisa ter mais 1 pq pode ser diferente para o Q predando outro host?
#taxa_disp -> lista com a taxa de dispersão de cada pop (PQNH)
#indv_migrante -> taxa de dispersão máxima que pode sair de cada população (PQNH/ 0Q0H)
#tempo_final -> qual vai ser a ultima geracao
#lin -> numero de linhas
#col -> numero de colunas
#pop_iniciais -> uma lista contendo as pop iniciais de cada ponto q tem pop inicial [[P,Q,N,H],[P,Q,N,H]], sendo q
	#pop_iniciais[0] vai ser estar no patch_iniciais[0]
#patch_iniciais -> uma lista contendo os patchs que ira iniciar com alguma pop [[l,c],[l,c]]
#viz-> é a lista dos vizinhos e distancias criada com a funcao criando_lista_viz(). É a primeira lista criada pela funcao citada.
#lista_raios-> é a lista de raios q tem no grid escolhido, é a segunda lista craida pela funcao criando_lista_viz
#TS= "total available searching time" ([P,Q1,Q2])
#TH = handling time ([P,Q1,Q2]) (tem dois Q pois o Q pega 2 host, então pode ter tempo de manuseio diferente para ambos. O q1 é para o host compartilhado e o Q2 para o host exclusivo do generalista.)
#porcentagem_popinihost= a porcemtagem de patches q cada host vai ocupar no inicio. ([N,H])
#inicial_h = quanto é a pop inicial do H (3)
#inicial_host1 = quanto é ap op inicialdo N (2)
#Pnumerico = media de qnts cotesia adultos emergem de uma única pulpa de diatraea
#Q1numerico = media de qnts tetrastichus adultos emergem de uma única pulpa de diatraea
#Q2numerico = media de qnts tetrastichus adultos emergem de uma única pulpa de HOST ALTERNATIVO
 


def modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico):
	qnd_onde_liberou=[]
	ocupacao_total=[[],[],[],[]] #criando a matriz de ocupacao ocupacao_total[pop][t]
	ocupacao_total_tempo=[0,0,0,0] #ocupacao_total_tempo[pop]

	media_regional=[[],[],[],[]] #criando a lista que vai acomodar a media regional media_regional[populacao][t]
	media_regional_tempo=[0,0,0,0] #criando a lista que vai acomodar a media regional media_regional[populacao][t]

	grid=lin*col #criando o valor de grid para poder fazer a media regional

	pop_migrante_tempo=[0,0,0,0]
	media_regional_migracao=[[0],[0],[0],[0]]

	#criando os vetores dos dados
	vetores=criador_vetor(lin,col,tempo_final)
	g_all=vetores[0]
	pop_migrante=vetores[1]
	visitacao_patchs=vetores[2]
	imigracao_patchs=vetores[3]
	emigracao_patchs=vetores[4]

	###########
	#INICIO DINAMICA POPULACIONAL
	###########

	for t in range(tempo_final):
		if t>100 and ocupacao_total[0][t-1]==0 and ocupacao_total[1][t-1]==0:
			tempo_gall_para=[]
			tempo_gall_host1=[]
			tempo_gall_host2=[]
			tempo_ocupacao_para=[]
			tempo_ocupacao_host1=[]
			tempo_ocupacao_host2=[]
			tempo_migracao_host1=[] 
			tempo_migracao_host2=[] 
			tempo_migracao_para=[]
			tempo_imigracao_host1=[]
			tempo_imigracao_host2=[]
			tempo_imigracao_para=[]
			nova_pop_h1=0
			nova_pop_h2=0
			if ocupacao_total[2][t-1]>0:
				nova_pop_h1=10
			if ocupacao_total[3][t-1]>0:
				nova_pop_h2=10

			for t_novo in range(tempo_final-t): 
				tempo_gall_para.append(0)
				tempo_gall_host1.append(nova_pop_h1)
				tempo_gall_host2.append(nova_pop_h2)
				tempo_ocupacao_para.append(0)
				tempo_ocupacao_host1.append(nova_pop_h1)
				tempo_ocupacao_host2.append(nova_pop_h2)
				tempo_migracao_host1.append(nova_pop_h1)
				tempo_migracao_host2.append(nova_pop_h2)
				tempo_migracao_para.append(0)
				tempo_imigracao_host1.append(nova_pop_h1)
				tempo_imigracao_host2.append(nova_pop_h2)
				tempo_imigracao_para.append(0)
			for l in range(lin):
				for c in range(col):
					g_all[l][c][0]=g_all[l][c][0]+tempo_gall_para
					g_all[l][c][1]=g_all[l][c][1]+tempo_gall_para
					g_all[l][c][2]=g_all[l][c][2]+tempo_gall_host1
					g_all[l][c][3]=g_all[l][c][3]+tempo_gall_host2

					pop_migrante[l][c][0]=pop_migrante[l][c][0]+tempo_migracao_para
					pop_migrante[l][c][1]=pop_migrante[l][c][1]+tempo_migracao_para
					pop_migrante[l][c][2]=pop_migrante[l][c][2]+tempo_migracao_host1
					pop_migrante[l][c][3]=pop_migrante[l][c][3]+tempo_migracao_host2


					imigracao_patchs[l][c][0]=imigracao_patchs[l][c][0]+tempo_imigracao_para
					imigracao_patchs[l][c][1]=imigracao_patchs[l][c][1]+tempo_imigracao_para
					imigracao_patchs[l][c][2]=imigracao_patchs[l][c][2]+tempo_imigracao_host1
					imigracao_patchs[l][c][3]=imigracao_patchs[l][c][3]+tempo_imigracao_host2
			ocupacao_total[0]=ocupacao_total[0]+tempo_ocupacao_para
			ocupacao_total[1]=ocupacao_total[1]+tempo_ocupacao_para
			ocupacao_total[2]=ocupacao_total[2]+tempo_ocupacao_host1
			ocupacao_total[3]=ocupacao_total[3]+tempo_ocupacao_host2	
			salva_arquivos(zzzz,g_all,pop_migrante,ocupacao_total,media_regional,media_regional_migracao,visitacao_patchs,imigracao_patchs,emigracao_patchs,a1,a2,a3,lambda1,lambda2,taxa_disp,indv_migrante,tempo_final,lin,col,qnd_onde_liberou)
			break

		for l in range(lin):
			for c in range(col):
				
				if t>0:
					tempo_dinamica=dinamica_iterativa(g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante,c,l,t,a1,a2,a3,lambda1,lambda2,TS,TH,H0,h0,f0,ocupacao_total,Pnumerico,Q1numerico)
					#############
					#migração
					#############
					tempo2=criando_popmigrante_densidade_patch(l,c,t,H0,h0,f0,indv_migrante,g_all,pop_migrante,emigracao_patchs)#criei a pop migrante do tempo
					pop_migrante=tempo2[0]
					emigracao_patchs=tempo2[1]
					for i,j in enumerate(pop_migrante[l][c]): #estou diminuindo a pop do patch pela pop q vai migrar
						g_all[l][c][i][t]-=j[t]
				else:
					tempo_dinamica=dinamica_tempo0(g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante,c,l,t,patch_iniciais)
				g_all=tempo_dinamica[0]
				pop_iniciais=tempo_dinamica[1]
				ocupacao_total_tempo=tempo_dinamica[2]
				media_regional_tempo=tempo_dinamica[3]
				pop_migrante=tempo_dinamica[4]

		#saiu do loop pela matriz, agora, no mesmo tempo, vem o segundo loop pela matriz para poder migrar:
		if t!=0: #isso eh para nao ter migracao no tempo zero
			#loop de lin e col
			for l in range(lin):
				for c in range(col):

					tempo=sorteio_vizinhos(g_all,pop_migrante,taxa_disp,viz,l,c,t,visitacao_patchs,imigracao_patchs,lista_distancia)
					g_all=tempo[0]
					visitacao_patchs=tempo[1]
					imigracao_patchs=tempo[2]
					pop_migrante_tempo=[0,0,0,0]
					for i,j in enumerate(pop_migrante[l][c]): #i é posição, j elemento
						pop_migrante_tempo[i]+=pop_migrante[l][c][i][t]

			for i,j in enumerate(pop_migrante[l][c]): #pondo a media regional
				media_regional_migracao[i].append(pop_migrante_tempo[i]/grid)

		###Bloco de checar se as pop possuem recurso
		if t!=0: #para soh acontecer dps do tempo 0
			for l in range(lin):
				for c in range(col):
					g_all=conferir_recurso_parasitoides(g_all,l,c,t)
		
		

		
		
		if t!=0: #aqui atualiza a media regional e ocupacao
			for l in range(lin):
				for c in range(col):
					if g_all[l][c][0][t]>0:
						ocupacao_total_tempo[0]+=1
						media_regional_tempo[0]+=g_all[l][c][0][t]
					if g_all[l][c][1][t]>0:
						ocupacao_total_tempo[1]+=1
						media_regional_tempo[1]+=g_all[l][c][1][t]
					if g_all[l][c][2][t]>0:
						ocupacao_total_tempo[2]+=1
						media_regional_tempo[2]+=g_all[l][c][2][t]
					if g_all[l][c][3][t]>0:
						ocupacao_total_tempo[3]+=1
						media_regional_tempo[3]+=g_all[l][c][3][t]


		ocupacao_total[0].append(ocupacao_total_tempo[0])
		ocupacao_total[1].append(ocupacao_total_tempo[1])
		ocupacao_total[2].append(ocupacao_total_tempo[2])
		ocupacao_total[3].append(ocupacao_total_tempo[3])

		media_regional[0].append(media_regional_tempo[0]/grid)
		media_regional[1].append(media_regional_tempo[1]/grid)
		media_regional[2].append(media_regional_tempo[2]/grid)
		media_regional[3].append(media_regional_tempo[3]/grid)


		media_regional_tempo=[0,0,0,0]
		ocupacao_total_tempo=[0,0,0,0]
			
	salva_arquivos(zzzz,g_all,pop_migrante,ocupacao_total,media_regional,media_regional_migracao,visitacao_patchs,imigracao_patchs,emigracao_patchs,a1,a2,a3,lambda1,lambda2,taxa_disp,indv_migrante,tempo_final,lin,col,qnd_onde_liberou)
	return
def criando_lista_patches_cada_hectare():
	lin=50
	col=50
	hectare1=[]
	hectare2=[]
	hectare3=[]
	hectare4=[]
	hectare5=[]
	hectare6=[]
	hectare7=[]
	hectare8=[]
	hectare9=[]
	hectare10=[]
	for l in range(lin):
		for c in range(col):
			entrou_em_qnts=0
			if l<10 and c<25:
				hectare1.append([l,c])
				entrou_em_qnts+=1
			elif l>9 and l<20 and c<25:
				hectare2.append([l,c])
				entrou_em_qnts+=1
			elif l>19 and l<30 and c<25:
				hectare3.append([l,c])
				entrou_em_qnts+=1
			elif l>29 and l<40 and c<25:
				hectare4.append([l,c])
				entrou_em_qnts+=1
			elif l>39 and l<50 and c<25:
				hectare5.append([l,c])
				entrou_em_qnts+=1
			elif l<10 and c>24:
				hectare6.append([l,c])
				entrou_em_qnts+=1
			elif l>9 and l<20 and c>24:
				hectare7.append([l,c])
				entrou_em_qnts+=1
			elif l>19 and l<30 and c>24:
				hectare8.append([l,c])
				entrou_em_qnts+=1
			elif l>29 and l<40 and c>24:
				hectare9.append([l,c])
				entrou_em_qnts+=1
			elif l>39 and l<50 and c>24:
				hectare10.append([l,c])
				entrou_em_qnts+=1
	arquivoviz=f'50x50_patches_por_hectare.txt'
	a=open(arquivoviz,"ab")		 
	dados_em_pck=pck.dumps([hectare1,hectare2,hectare3,hectare4,hectare5,hectare6,hectare7,hectare8,hectare9,hectare10])   
	a.write(dados_em_pck)
	a.close()
	return [hectare1,hectare2,hectare3,hectare4,hectare5,hectare6,hectare7,hectare8,hectare9,hectare10]
def criando_lista_patches_cada_borda():
	lin=50
	col=50
	hectare1=[]
	hectare2=[]
	hectare3=[]
	hectare4=[]
	hectare5=[]
	hectare6=[]
	hectare7=[]
	hectare8=[]
	hectare9=[]
	hectare10=[]
	for l in range(lin):
		for c in range(col):
			entrou_em_qnts=0
			if (l<10 and c==24) or (l<10 and c==0) or (l==0 and c<25) or (l==9 and c<25):
				hectare1.append([l,c])
				entrou_em_qnts+=1
			elif (l>9 and l<20 and c==24) or (l>9 and l<20 and c==0) or (l==10 and c<25) or (l==19 and c<25):
				hectare2.append([l,c])
				entrou_em_qnts+=1
			elif (l>19 and l<30 and c==24) or (l>19 and l<30 and c==0) or (l==20 and c<25) or (l==29 and c<25):
				hectare3.append([l,c])
				entrou_em_qnts+=1
			elif (l>29 and l<40 and c==24) or (l>29 and l<40 and c==0) or (l==30 and c<25) or (l==39 and c<25):
				hectare4.append([l,c])
				entrou_em_qnts+=1
			elif (l>39 and l<50 and c==24) or (l>39 and l<50 and c==0) or (l==40 and c<25) or (l==49 and c<25):
				hectare5.append([l,c])
				entrou_em_qnts+=1
			elif (l<10 and c==25) or (l<10 and c==49) or (l==0 and c>24) or (l==9 and c>24):
				hectare6.append([l,c])
				entrou_em_qnts+=1
			elif (l>9 and l<20 and c==25) or (l>9 and l<20 and c==49) or (l==10 and c>24) or (l==19 and c>24):
				hectare7.append([l,c])
				entrou_em_qnts+=1
			elif (l>19 and l<30 and c==25) or (l>19 and l<30 and c==49) or (l==20 and c>24) or (l==29 and c>24):
				hectare8.append([l,c])
				entrou_em_qnts+=1
			elif (l>29 and l<40 and c==25) or (l>29 and l<40 and c==49) or (l==30 and c>24) or (l==39 and c>24):
				hectare9.append([l,c])
				entrou_em_qnts+=1
			elif (l>39 and l<50 and c==25) or (l>39 and l<50 and c==49) or (l==40 and c>24) or (l==49 and c>24):
				hectare10.append([l,c])
				entrou_em_qnts+=1
	arquivoviz=f'50x50_patches_por_hectare_bordas.txt'
	a=open(arquivoviz,"ab")		 
	dados_em_pck=pck.dumps([hectare1,hectare2,hectare3,hectare4,hectare5,hectare6,hectare7,hectare8,hectare9,hectare10])   
	a.write(dados_em_pck)
	a.close()
	return [hectare1,hectare2,hectare3,hectare4,hectare5,hectare6,hectare7,hectare8,hectare9,hectare10]
def criando_lista_distancia(lin,col,distancia_maxima):
	lista_distancia=[]
	for l in range(lin):
		for c in range(col):
			if l==0 and c==0:
				continue
			else:
				distancia=(((0-l)**2)+((0-c)**2))**(1/2)
				if distancia<=distancia_maxima:
					lista_distancia.append(distancia)
	lista_distancia=sorted(set(lista_distancia)) #qnd o loop termina, eu retiro os valores repetidos e ordeno a lista
	return lista_distancia
def checando_existencia_lista_distancias(lin,col,distancia_maxima):
	arquivodist=f'{lin}x{col}_max{distancia_maxima}_lista_distancia.txt'
	try:
		distancia_file=open(arquivodist,"rb")
		distancia_pck=distancia_file.read()
		distancia_file.close()
		lista_distancia=pck.loads(distancia_pck)
	except FileNotFoundError:
		print("Não achou a lista de distância. Iniciando processo de criação.")
		lista_distancia=criando_lista_distancia(lin,col,distancia_maxima)
		a=open(arquivodist,"ab")		 
		dados_em_pck=pck.dumps(lista_distancia)   
		a.write(dados_em_pck)
		a.close()
	return lista_distancia
def checando_existencia_lista_viz_e_lista_de_patches(lin,col,distancia_maxima):
	arquivoviz=f'{lin}x{col}vizinhos.txt'
	try:
		viz_file=open(arquivoviz,"rb")
		viz_pck=viz_file.read()
		viz_file.close()
		viz=pck.loads(viz_pck)
	except FileNotFoundError:
		print("Não achou a lista de vizinhos. Iniciando processo de criação.")
		viz=criando_lista_viz(lin,col,distancia_maxima)
		a=open(arquivoviz,"ab")		 
		dados_em_pck=pck.dumps(viz)   
		a.write(dados_em_pck)
		a.close()
	######aqui vai checar borda
	arquivoviz=f'50x50_patches_por_hectare_bordas.txt'
	try:
		viz_file=open(arquivoviz,"rb")
		viz_pck=viz_file.read()
		viz_file.close()
		lista_borda=pck.loads(viz_pck)
	except FileNotFoundError:
		print("Não achou a lista de patches na borda. Iniciando processo de criação.")
		lista_borda=criando_lista_patches_cada_borda()
		a=open(arquivoviz,"ab")		 
		dados_em_pck=pck.dumps(lista_borda)   
		a.write(dados_em_pck)
		a.close()
	#aqui vai checar lista total
	arquivoviz=f'50x50_patches_por_hectare.txt'
	try:
		viz_file=open(arquivoviz,"rb")
		viz_pck=viz_file.read()
		viz_file.close()
		lista_patches=pck.loads(viz_pck)
	except FileNotFoundError:
		print("Não achou a lista de patches do hectare. Iniciando processo de criação.")
		lista_patches=criando_lista_patches_cada_hectare()
		a=open(arquivoviz,"ab")		 
		dados_em_pck=pck.dumps(lista_patches)   
		a.write(dados_em_pck)
		a.close()
	return viz,lista_patches,lista_borda
def criando_lista_viz(lin,col,distancia_maxima):
	lista_raios=[]
	viz=[] #criando uma matriz vazia
	###bloco de criar a matriz de vizinhos vazia para poder adicionar os vizinhos
	#e como ja vai percorrer uma vez a matriz inteira, eu aproveito e ja calculo todos os vizinhos possiveis e faco a lista
	#de vizinhos possiveis.
	for l in range(lin): #percorrendo os loops de linhas
		viz.append([]) #colocando as linhas na matriz
		for c in range(col): #percorrendo os loops de colunas
			viz[l].append([]) #colocando as colunas na matriz
			if l==0 and c==0: #aqui checo soh pra ver se eh o primeiro patch de todos
				continue #se for ele nao faz a conta, pois ele eh raio zero
			else:	#qnd n for o patch 0,0 ele vai entrar aqui
				distancia_tempo=((0-l)**2+(0-c)**2)**(1/2)
				if distancia_tempo<=distancia_maxima:
					lista_raios.append(distancia_tempo)

					 #aqui eu faco a conta da distancia e adiciona na lista
	lista_raios=sorted(set(lista_raios)) #qnd o loop termina, eu retiro os valores repetidos e ordeno a lista
	#pois agora ela vai ter todas as distancias presentes nesse grid, e ordenadas por raio.

	###bloco de adicionar listas dos raios possiveis dentro da matriz
	for l in range(lin):
		for c in range(col):
			for i in range(len(lista_raios)):
				viz[l][c].append([]) #colocando todos os raios possiveis da matriz

	###bloco de achar os vizinhos e por no raio correspondente
	for l in range(lin): #loop da linha para o patch focal
		for c in range(col): #loop da coluna para o patch focal
			for ll in range(-lin,lin): #loop de linha do vizinho do patch focal
				for cc in range(-col,col): #loop de coluna do vizinho do patch focal
					if l==ll and c==cc: #checando se o patch focal eh o mesmo q o vizinho
						continue #se for, acontece nada
					else: 
						dist=((l-ll)**2+(c-cc)**2)**(1/2) #calculando a distancia do patch focal(l,c) pro vizinho(ll,cc)
						for i,j in enumerate(lista_raios): #loop para rodar nos raios existentes, pegando a casa q ele esta(i) e o valor do raio (j)
							if dist==j: #checando se o valor da distancia calculado eh igual a algum raio e descobrindo em qual lista de vizinho por
								if ll>=lin:
									tempolin=(lin-1)+((lin-1)-ll)
								elif ll<0:
									tempolin=ll*(-1)
								else:
									tempolin=ll
								if cc>=col:
									tempocol=(col-1)+((col-1)-cc)
								elif cc<0:
									tempocol=cc*(-1)
								else:
									tempocol=cc
								viz[l][c][i].append([tempolin,tempocol])
								break										
	return viz,lista_raios
def criando_popmigrante_densidade_patch(l,c,t,H0,h0,f0,indv_migrante,g_all,pop_migrante,emigracao_patchs):
	#p
	if g_all[l][c][2][t]==0:
		pop_migrante[l][c][0].append(int(g_all[l][c][0][t]))
		emigracao_patchs[l][c][0].append(int(g_all[l][c][0][t])) 
	else:
		pop_migrante[l][c][0].append(int(indv_migrante[0]*(H0[0]/(H0[0]+g_all[l][c][2][t]))*((g_all[l][c][0][t]**2)/(g_all[l][c][0][t]+f0[0]))))
		emigracao_patchs[l][c][0].append(int(indv_migrante[0]*(H0[0]/(H0[0]+g_all[l][c][2][t]))*((g_all[l][c][0][t]**2)/(g_all[l][c][0][t]+f0[0]))))
	#q
	if g_all[l][c][2][t]==0 and g_all[l][c][3][t]==0:
		pop_migrante[l][c][1].append(int(g_all[l][c][1][t]))
		emigracao_patchs[l][c][1].append(int(g_all[l][c][1][t]))
	else:
		pop_migrante[l][c][1].append(int(indv_migrante[1]*(H0[1]/(H0[1]+(g_all[l][c][2][t]+g_all[l][c][3][t])))*((g_all[l][c][1][t]**2)/(g_all[l][c][1][t]+f0[1]))))
		emigracao_patchs[l][c][1].append(int(indv_migrante[1]*(H0[1]/(H0[1]+(g_all[l][c][2][t]+g_all[l][c][3][t])))*((g_all[l][c][1][t]**2)/(g_all[l][c][1][t]+f0[1]))))
	#n
	pop_migrante[l][c][2].append(int((indv_migrante[2]*(g_all[l][c][2][t]**2))/(g_all[l][c][2][t]+h0[0])))
	emigracao_patchs[l][c][2].append(int((indv_migrante[2]*(g_all[l][c][2][t]**2))/(g_all[l][c][2][t]+h0[0])))
	#h 
	pop_migrante[l][c][3].append(int((indv_migrante[3]*(g_all[l][c][3][t]**2))/(g_all[l][c][3][t]+h0[1])))
	emigracao_patchs[l][c][3].append(int((indv_migrante[3]*(g_all[l][c][3][t]**2))/(g_all[l][c][3][t]+h0[1])))

	return pop_migrante,emigracao_patchs
def sorteio_vizinhos(g_all,pop_migrante,taxa_disp,viz,l,c,t,visitacao_patchs,imigracao_patchs,lista_distancia):
	#criando a pop migrante temporaria para poder diminuir desta variavel ateh chegar a zero
	pop_migrante_tempo=[]
	for i,p in enumerate(taxa_disp):
		pop_migrante_tempo.append(pop_migrante[l][c][i][t])
	for i,p in enumerate(pop_migrante_tempo): #i eh a casa. p eh o elemento #loop para rodar a qnt certa de pop e ja usar as taxa dela de dispersao
		x=0
		contador_para=0

		while pop_migrante_tempo[i]>0 and x<len(lista_distancia) and contador_para<3:
			if i==0 or i==1:
				contador_para+=1
			for j in viz[l][c][x]:
				if (i==0 or i==1) and x>=2:
					break
				sorteado=random.choice(viz[l][c][x])

				indv=(taxa_disp[i]/lista_distancia[x])*pop_migrante[l][c][i][t]
				if indv<1: #conferindo se vai migrar um valor menor do que 1, se for eu arredondo para 1
					indv=1
				if indv>pop_migrante_tempo[i]:
					indv=pop_migrante_tempo[i]
				pop_migrante_tempo[i]-=indv
				indv=int(indv)
				g_all[sorteado[0]][sorteado[1]][i][t]+=indv
				visitacao_patchs[sorteado[0]][sorteado[1]][i][t]+=1
				imigracao_patchs[sorteado[0]][sorteado[1]][i][t]+=indv
				if pop_migrante_tempo[i]<=0:
					break
			x+=1
	return g_all,visitacao_patchs,imigracao_patchs
def dinamica_tempo0(g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante,c,l,t,patch_iniciais):
	achou_patch=0
	for i,j in enumerate(patch_iniciais): #aqui eu fiz um for para passar pela lista de patchs iniciais para
		if j[0]==l and j[1]==c: #eu poder comparar com qual linha/coluna está. se alguma linha/coluna for igual ao da lista
			#de patchs iniciais, ele adiciona a pop inicial. o 'i' serve para saber qual casa o patch esta na lista patch_iniciais
			#e ai a mesma casa é usada para achar a pop inicial daquele patch.
			
			#pop 0 - P
			achou_patch=1
			g_all[l][c][0].append(pop_iniciais[i][0])
			if g_all[l][c][0][t]>0.01:
				ocupacao_total_tempo[0]+=1
				media_regional_tempo[0]+=g_all[l][c][0][t]


			#pop 1 - Q
			#testando
			g_all[l][c][1].append(pop_iniciais[i][1])
			if g_all[l][c][1][t]>0.01:
				ocupacao_total_tempo[1]+=1
				media_regional_tempo[1]+=g_all[l][c][1][t]


			#pop 2 - N
			g_all[l][c][2].append(pop_iniciais[i][2])
			if g_all[l][c][2][t]>0.01:
				ocupacao_total_tempo[2]+=1
				media_regional_tempo[2]+=g_all[l][c][2][t]
			
			
			#pop 3 - H
			g_all[l][c][3].append(pop_iniciais[i][3])
			if g_all[l][c][3][t]>0.01:
				ocupacao_total_tempo[3]+=1
				media_regional_tempo[3]+=g_all[l][c][3][t]


			#colocando zero na migracao
			pop_migrante[l][c][0].append(0)
			pop_migrante[l][c][1].append(0)
			pop_migrante[l][c][2].append(0)
			pop_migrante[l][c][3].append(0)

	if achou_patch==0:
		#colocando zero nas pop iniciais
		g_all[l][c][0].append(0)
		g_all[l][c][1].append(0)
		g_all[l][c][2].append(0)
		g_all[l][c][3].append(0)

		#colocando zero na migracao
		pop_migrante[l][c][0].append(0)
		pop_migrante[l][c][1].append(0)
		pop_migrante[l][c][2].append(0)
		pop_migrante[l][c][3].append(0)
	return g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante
def dinamica_iterativa(g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante,c,l,t,a1,a2,a3,lambda1,lambda2,TS,TH,H0,h0,f0,ocupacao_total,Pnumerico,Q1numerico):
	
	rfP=(math.exp((-a1*TS[0]*g_all[l][c][0][t-1])/(1+a1*TH[0]*g_all[l][c][2][t-1])))#resposta funcional especialista(P) (P para o N)
	rfQ1=(math.exp((-a2*TS[1]*g_all[l][c][1][t-1])/(1+a2*TH[1]*g_all[l][c][2][t-1])))#resposta funcional generalsita(Q) (Q para o N)
	rfQ2=(math.exp((-a3*TS[2]*g_all[l][c][1][t-1])/(1+a3*TH[2]*g_all[l][c][3][t-1])))#resposta funcional generalsita(Q) (Q para o H) #to-do: conferir a resposta funcional do generalista com o host alternativo (q->h)
	
	#vou checar se a resposta funcional é zero ou menor, pq se for a equação dos bixos daria zero pq multiplica ela
	if rfP<=0:
		rfP=1 #coloco igual a 1 pq ai a resposta funcional nao vai alterar em nada as equações.
	if rfQ1<=0:
		rfQ1=1
	if rfQ2<=0:
		rfQ2=1

	#bloquinho dinamica populacao 0 (P)
	if g_all[l][c][0][t-1]==0:
		g_all[l][c][0].append(0)
	else:		
		g_all[l][c][0].append(g_all[l][c][2][t-1]*(1-rfP)*Pnumerico)
		if g_all[l][c][0][t]<=0.01:
			g_all[l][c][0][t]=0

	#bloquinho dinamica populacao 1 (Q):
	#teste
	if g_all[l][c][1][t-1]==0: 
		g_all[l][c][1].append(0)
	else:
		g_all[l][c][1].append(g_all[l][c][2][t-1]*rfP*(1-rfQ1)*Q1numerico+g_all[l][c][3][t-1]*(1-rfQ2)*Q1numerico) #to-do: essa equação só leva em conta 1 hospedeiro.
		if g_all[l][c][1][t]<=0.01:
			g_all[l][c][1][t]=0

	#bloquinho dinamica populacao 2 (N)
	if g_all[l][c][2][t-1]==0:
		g_all[l][c][2].append(0)
	else:
		g_all[l][c][2].append(g_all[l][c][2][t-1]*lambda1*rfP*rfQ1)
		if g_all[l][c][2][t]<=0.01:
			g_all[l][c][2][t]=0
		if g_all[l][c][2][t]>1000 and ocupacao_total[0][t-1]==0 and ocupacao_total[1][t-1]==0: 
			g_all[l][c][2][t]=1000


	#bloquinho dinamica populacao 3 (h):
	#teste:
	if g_all[l][c][3][t-1]==0:
		g_all[l][c][3].append(0)
	else:
		g_all[l][c][3].append(g_all[l][c][3][t-1]*lambda2*rfQ2) 
		if g_all[l][c][3][t]<=0.01:
			g_all[l][c][3][t]=0
		if g_all[l][c][3][t]>1000 and ocupacao_total[1][t-1]==0:
			g_all[l][c][3][t]=1000
	return g_all,pop_iniciais,ocupacao_total_tempo,media_regional_tempo,pop_migrante
def conferir_recurso_parasitoides(g_all,l,c,t):
	if g_all[l][c][0][t]>0 and g_all[l][c][2][t]==0:
		g_all[l][c][0][t]=0
	if g_all[l][c][1][t]>0 and g_all[l][c][2][t]==0 and g_all[l][c][3][t]==0:
		g_all[l][c][1][t]=0
	return g_all
def salva_arquivos(zzzz,g_all,pop_migrante,ocupacao_total,media_regional,media_regional_migracao,visitacao_patchs,imigracao_patchs,emigracao_patchs,a1,a2,a3,lambda1,lambda2,taxa_disp,indv_migrante,tempo_final,lin,col,qnd_onde_liberou):
	

	nome=f'{zzzz}_qnd_onde_liberou.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(qnd_onde_liberou)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_g_all.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(g_all)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_pop_migrante.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(pop_migrante)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_ocupacao_total.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(ocupacao_total)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_media_regional.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(media_regional)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_media_regional_migracao.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(media_regional_migracao)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_visitacao_patchs.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(visitacao_patchs)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_imigracao_patchs.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(imigracao_patchs)   
	a.write(dados_em_pck)
	a.close()

	nome=f'{zzzz}_emigracao_patchs.txt'
	a=open(nome,"ab")		 
	dados_em_pck=pck.dumps(emigracao_patchs)   
	a.write(dados_em_pck)
	a.close()

	nome=f'read_me.txt'
	a=open(nome,"w")
	conteudo_readme=f'A ordem do nome é tipo.do.arquivo.txt\n \n \n Os arquivos possuem listas:\n g_all[linha][coluna][populacao][tempo]: é a populacao de cada patch por tempo\npop_migrante[linha][coluna][populacao][tempo]: é o vetor de migracao\nocupacao_total[populacao][tempo]: é a ocupacao de patchs por tempo\nmedia_regional[populacao][tempo]: é a media regional de populacao \nmedia_regional_migracao[populacao][tempo]: é a media regional de migracao\nvisitacao_patchs[linha][coluna][populacao][tempo]:é quais patchs receberam imigrantes\nimigracao_patchs[linha][coluna][populacao][tempo]: é o quanto cada patch recebeu de migrantes, ou seja, a imigracao\nemigracao_patchs[linha][coluna][pop][tempo]:é o qnt que cada pop de cada patch contribuiu para a migraçao, ou seja, a emigraçao \n \nOs valores desse conjunto de simulação são:\na1={a1}\na2={a2}\na3={a3}\nlambda1={lambda1}\nlambda2={lambda2}\ntaxa_disp={taxa_disp}\nindv_migrante={indv_migrante}\ntempo_final={tempo_final}\nlin={lin}\ncol={col}'
	a.write(conteudo_readme)
	a.close()
	return
def criador_vetor(lin,col,tempo_final):
	#Primeiro eu crio uma lista vazia, e dps percorro o numero de linhas e vou adicionando uma lista vazia opr iteracao.
	#Em seguida, dentro de cada linha, eu percorro o numero de colunas e vou adicionando lista vazia por coluna, por linha
	#e em seguida percorro o numero de pop dentro de cada coluna e adiciona listas vazias iguais o numero de pop
	#dentro de cada coluna dentro de cada linha.
	g_all=[] #g_all[lin][col][pop][tempo]
	pop_migrante=[] #pop_migrante[lin][col][pop][tempo]
	visitacao_patchs=[]
	imigracao_patchs=[]
	emigracao_patchs=[]
	for l in range(lin):
		g_all.append([])
		pop_migrante.append([])
		visitacao_patchs.append([])
		imigracao_patchs.append([])
		emigracao_patchs.append([])
		for c in range(col):
			g_all[l].append([[],[],[],[]])
			pop_migrante[l].append([[],[],[],[]])
			visitacao_patchs[l].append([[],[],[],[]])
			imigracao_patchs[l].append([[],[],[],[]])
			emigracao_patchs[l].append([[0],[0],[0],[0]])
			for t in range(tempo_final):
				visitacao_patchs[l][c][0].append(0)
				visitacao_patchs[l][c][1].append(0)
				visitacao_patchs[l][c][2].append(0)
				visitacao_patchs[l][c][3].append(0)

				imigracao_patchs[l][c][0].append(0)
				imigracao_patchs[l][c][1].append(0)
				imigracao_patchs[l][c][2].append(0)
				imigracao_patchs[l][c][3].append(0)
	return g_all,pop_migrante,visitacao_patchs,imigracao_patchs,emigracao_patchs
def sorteio_pop_inicial_host(lin,col,porcentagem_popinihost,inicial_host1,inicial_host2,solt_valor,loop,lista_patches): #retorna em primeiro os patches iniciais e em segundo as populações iniciais APENAS DO HOST
	pop_iniciais=[]
	grid=lin*col
	#% de n:
	numero_patch_lib_n=grid*porcentagem_popinihost[0]
	numero_patch_lib_n/=10
	#% de h:
	numero_patch_lib_h=grid*porcentagem_popinihost[1]
	numero_patch_lib_h/=10

	if loop==25 or loop==28: #to-do: mudar isso para quando fizer os loops
		porcentagem_n=int(numero_patch_lib_n-(solt_valor*2)) #n
		porcentagem_h=int(numero_patch_lib_h-(solt_valor*2)) #h #aqui ta indo com os valores do outro host
	else:
		porcentagem_n=int(numero_patch_lib_n-solt_valor) #n #to-do: mudar isso para quando fizer os loops
		porcentagem_h=int(numero_patch_lib_h-solt_valor) #h #aqui ta indo com os valores do outro host
	patch_ini=[]

	#pondo n:
	for posicao_patch,lista_dos_patches in enumerate(lista_patches):
		for i in range(porcentagem_n):
			patch_tempo=random.choice(lista_dos_patches)
			while patch_tempo in patch_ini:
				patch_tempo=random.choice(lista_dos_patches)
			patch_ini.append(patch_tempo)
			pop_iniciais.append([0,0,inicial_host1,0])
	#pondo h:
	for posicao_patch,lista_dos_patches in enumerate(lista_patches):
		for i in range(porcentagem_h):
			patch_tempo=random.choice(lista_dos_patches)

			while patch_tempo in patch_ini and pop_iniciais[patch_ini.index(patch_tempo)][3]>0:
				patch_tempo=random.choice(lista_dos_patches)
			if patch_tempo in patch_ini:
				pop_iniciais[patch_ini.index(patch_tempo)][3]=inicial_host2
			else:
				patch_ini.append(patch_tempo)
				pop_iniciais.append([0,0,0,inicial_host2])
	return patch_ini,pop_iniciais
def alocando_patch_central_com_todas_pop(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen): #retorna em primeiro os patches iniciais e em segundo as populações iniciais APENAS DO HOST
	pop_iniciais=[]
	patch_ini=[]
	grid=lin*col
	meio=[int(lin/2-1),int(col/2-1)]

	pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,inicial_host2])
	patch_ini.append(meio)
	return patch_ini,pop_iniciais
def alocar_parasitoide(borda_meio,solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,hec10_bordas,inicial_host1,inicial_host2,para_lista_inicial,loop):
	if borda_meio==0:
		patches_sorteados=[[],[]]
		tempo=BORDA_alocando_parasitoide(solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,lista_borda,inicial_host1,inicial_host2,para_lista_inicial,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		return patch_iniciais,pop_iniciais,patches_sorteados

	else:
		tempo=MEIO_alocando_parasitoide(solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,inicial_host1,inicial_host2,para_lista_inicial,loop,lista_patches)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		patches_sorteados=tempo[2]
		return patch_iniciais,pop_iniciais,patches_sorteados
def MEIO_alocando_parasitoide(solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,inicial_host1,inicial_host2,para_lista_inicial,loop,lista_patches):
	if loop==23 or loop==24 or loop==26 or loop==27: #to-do: arrumar isso quando fizer os loops
		patches_sorteados=[[],[]]
		for qual_hectare in lista_patches:
			for sorteio in range(solt_valor):
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],para_lista_inicial[1][solt_posi],inicial_host1,0])
				patches_sorteados[0].append(patch_tempo)
	elif loop==25 or loop==28:
		patches_sorteados=[[],[]]
		#pondo cotesia:
		for qual_hectare in lista_patches:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],0,inicial_host1,0])
				patches_sorteados[0].append(patch_tempo)
		#pondo tetrastichus
		for qual_hectare in lista_patches:
			for sorteio in range(solt_valor):
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([0,para_lista_inicial[1][solt_posi],inicial_host1,0])
				patches_sorteados[1].append(patch_tempo)
	
	
	return patch_iniciais,pop_iniciais,patches_sorteados
def BORDA_alocando_parasitoide(solt_posi,solt_valor,patch_iniciais,pop_iniciais,lin,col,lista_borda,inicial_host1,inicial_host2,para_lista_inicial,loop):
	
	if loop==19 or loop==20 or loop==23 or loop==24:
		patches_sorteados=[[],[]]
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],para_lista_inicial[1][solt_posi],inicial_host1,0])
				patches_sorteados[0].append(patch_tempo)
	elif loop==21 or loop==25:
		patches_sorteados=[[],[]]
		#pondo cotesia:
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],0,inicial_host1,0])
				patches_sorteados[0].append(patch_tempo)
		#pondo tetrastichus
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([0,para_lista_inicial[1][solt_posi],inicial_host1,0])
				patches_sorteados[1].append(patch_tempo)
	elif loop==22 or loop==26 or loop==27:
		patches_sorteados=[[],[]]
		#pondo cotesia:
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],0,inicial_host1,0])
				patches_sorteados[0].append(patch_tempo)
		#pondo segunda cotesia
		for qual_hectare in lista_borda:
			for sorteio in range(solt_valor):
				
				patch_tempo=random.choice(qual_hectare)
				while patch_tempo in patch_iniciais:
					
					patch_tempo=random.choice(qual_hectare)
				patch_iniciais.append(patch_tempo)
				pop_iniciais.append([para_lista_inicial[0][solt_posi],0,inicial_host1,0])
				patches_sorteados[0].append(patch_tempo)
	
	return patch_iniciais,pop_iniciais,patches_sorteados
def checando_necessidade_liberacao(t, g_all,lin,col,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,loop,patches_sorteados,borda_meio,lista_patches,lista_borda,qnd_onde_liberou):
	#lista_patches,lista_borda
	teve_liberacao=False
	pop_n_inicial=0
	pop_n_t=0
	for posicao_hectare,hectare in enumerate(lista_patches):
		pop_n_inicial=0
		pop_n_t=0
		for patch in hectare:
			pop_n_t+=g_all[patch[0]][patch[1]][2][t]
		if pop_n_t>=925: #se mudar a pop inicial tem q mudar isso.
			teve_liberacao=True
			qnd_onde_liberou.append(f't{t}h{posicao_hectare}')
			if borda_meio==0: #borda
				patches_para_liberar=lista_borda[posicao_hectare]
				g_all=nova_liberacao(t, g_all,lin,col,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,loop,patches_sorteados,borda_meio,patches_para_liberar)
			else: #meio
				patches_para_liberar=lista_patches[posicao_hectare]
				g_all=nova_liberacao(t, g_all,lin,col,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,loop,patches_sorteados,borda_meio,patches_para_liberar)
	return g_all,qnd_onde_liberou,teve_liberacao
def nova_liberacao(t, g_all,lin,col,solt_posi,solt_valor,hec10_bordas,para_lista_inicial,loop,patches_sorteados,borda_meio,patches_para_liberar):
	#aqui eu preciso garantir q qnd é separado eu GARANTA q nao acontece de ficar igual
	if loop==21 or loop==19:
		patches_q_um_para_foi=[]
		#pondo cotesia
		for sorteio in range(solt_valor):
			patch_tempo=random.choice(patches_para_liberar)
			while patch_tempo in patches_q_um_para_foi:
				patch_tempo=random.choice(patches_para_liberar)
			g_all[patch_tempo[0]][patch_tempo[1]][0][t]+=para_lista_inicial[0][solt_posi]
			patches_q_um_para_foi.append(patch_tempo)
		#pondo tetrastichus:
		for sorteio in range(solt_valor):
			patch_tempo=random.choice(patches_para_liberar)
			while patch_tempo in patches_q_um_para_foi:
				patch_tempo=random.choice(patches_para_liberar)
			patch_tempo=random.choice(patches_para_liberar)
			g_all[patch_tempo[0]][patch_tempo[1]][1][t]+=para_lista_inicial[1][solt_posi]
	elif loop==20 : 
		patches_q_um_para_foi=[]
		#pondo ambos
		for sorteio in range(solt_valor):
			patch_tempo=random.choice(patches_para_liberar)
			while patch_tempo in patches_q_um_para_foi:
				patch_tempo=random.choice(patches_para_liberar)
			g_all[patch_tempo[0]][patch_tempo[1]][0][t]+=para_lista_inicial[0][solt_posi]
			g_all[patch_tempo[0]][patch_tempo[1]][1][t]+=para_lista_inicial[1][solt_posi]
			patches_q_um_para_foi.append(patch_tempo)
	elif loop==22:
		patches_q_um_para_foi=[]
		#pondo primeira cotesia
		for sorteio in range(solt_valor):
			patch_tempo=random.choice(patches_para_liberar)
			while patch_tempo in patches_q_um_para_foi:
				patch_tempo=random.choice(patches_para_liberar)
			g_all[patch_tempo[0]][patch_tempo[1]][0][t]+=para_lista_inicial[0][solt_posi]
			patches_q_um_para_foi.append(patch_tempo)
		#pondo segunda cotesia:
		for sorteio in range(solt_valor):
			patch_tempo=random.choice(patches_para_liberar)
			while patch_tempo in patches_q_um_para_foi:
				patch_tempo=random.choice(patches_para_liberar)
			patch_tempo=random.choice(patches_para_liberar)
			g_all[patch_tempo[0]][patch_tempo[1]][0][t]+=para_lista_inicial[0][solt_posi]

	return g_all
def condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop):
	meio=[int(lin/2-1),int(col/2-1)]

	if loop==29:#100h1_0h2_q
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					pop_iniciais.append([0,0,inicial_host1,0])
	elif loop==30: #50h1_50h2_q
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l>meio[0]:
						pop_iniciais.append([0,0,0,inicial_host2])
					else:
						pop_iniciais.append([0,0,inicial_host1,0])	
	elif loop==31: #75h1_25h2_q
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l>meio[0] and c>meio[1]:
						pop_iniciais.append([0,0,0,inicial_host2])
					else:
						pop_iniciais.append([0,0,inicial_host1,0])
	elif loop==32: #25h1_75h2_q
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l<=meio[0] and c<=meio[1]:
						pop_iniciais.append([0,0,inicial_host1,0])
					else:
						pop_iniciais.append([0,0,0,inicial_host2])
	elif loop==33: #100%_100%
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				patch_ini.append([l,c])
				if l==meio[0] and c==meio[1]:
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1/2,inicial_host2/2])
				else:
					pop_iniciais.append([0,0,inicial_host1/2,inicial_host2/2])
	elif loop==34: #75h1_25h2_b
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l in [0,1,2,3,49,48,47] or c in [0,1,2,49,48,47]: #4 linhas de cada lado e 3 linhas cima e baixo
						pop_iniciais.append([0,0,0,inicial_host2])
					else:
						pop_iniciais.append([0,0,inicial_host1,0])
	elif loop==35: #50h1_50h2_b
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l>=7 and l<=42 and c>=7 and c<42: #4 linhas de cada lado e 3 linhas cima e baixo
						pop_iniciais.append([0,0,inicial_host1,0])
					else:
						pop_iniciais.append([0,0,0,inicial_host2])
	elif loop==36: #25h1_75h2_b
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l>=12 and l<=36 and c>=12 and c<=36: #4 linhas de cada lado e 3 linhas cima e baixo
						pop_iniciais.append([0,0,inicial_host1,0])
					else:
						pop_iniciais.append([0,0,0,inicial_host2])
	elif loop==37: #50h1_50h2_q2
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if (l>meio[0] and c>meio[1]) or (l<=meio[0] and c<=meio[1]):
						pop_iniciais.append([0,0,inicial_host1,0])
					else:
						pop_iniciais.append([0,0,0,inicial_host2])

	#apartir daqui sao os sempre100%h1
	elif loop==38: #100h1_75h2q
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l<=meio[0] and c<=meio[1]:
						pop_iniciais.append([0,0,inicial_host1,0])
					else:
						pop_iniciais.append([0,0,inicial_host1/2,inicial_host2/2])
	elif loop==39: #100h1_50h2q
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l>meio[0]:
						pop_iniciais.append([0,0,inicial_host1/2,inicial_host2/2])
					else:
						pop_iniciais.append([0,0,inicial_host1,0])
	elif loop==40: #100h1_50h2q2
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if (l>meio[0] and c>meio[1]) or (l<=meio[0] and c<=meio[1]):
						pop_iniciais.append([0,0,inicial_host1,0])
					else:
						pop_iniciais.append([0,0,inicial_host1/2,inicial_host2/2])
	elif loop==41: #100h1_25h2q
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l>meio[0] and c>meio[1]:
						pop_iniciais.append([0,0,inicial_host1/2,inicial_host2/2])
					else:
						pop_iniciais.append([0,0,inicial_host1,0])
	elif loop==42: #100h1_25h2b
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l in [0,1,2,3,49,48,47] or c in [0,1,2,49,48,47]: #4 linhas de cada lado e 3 linhas cima e baixo
						pop_iniciais.append([0,0,inicial_host1/2,inicial_host2/2])
					else:
						pop_iniciais.append([0,0,inicial_host1,0])
	elif loop==43: #100h1_50h2b
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l>=7 and l<=42 and c>=7 and c<42: #4 linhas de cada lado e 3 linhas cima e baixo
						pop_iniciais.append([0,0,inicial_host1,0])
					else:
						pop_iniciais.append([0,0,inicial_host1/2,inicial_host2/2])
	elif loop==44: #100h1_75h2b
		patch_ini=[]
		pop_iniciais=[]
		for l in range(lin):
			for c in range(col):
				if l==meio[0] and c==meio[1]:
					patch_ini.append([l,c])
					pop_iniciais.append([inicial_esp,inicial_gen,inicial_host1,0])
				else:
					patch_ini.append([l,c])
					if l>=12 and l<=36 and c>=12 and c<=36: #4 linhas de cada lado e 3 linhas cima e baixo
						pop_iniciais.append([0,0,inicial_host1,0])
					else:
						pop_iniciais.append([0,0,inicial_host1/2,inicial_host2/2])
	
	
	return patch_ini,pop_iniciais
def loops_cenarios(loop,viz,lista_patches,lista_borda,inicial_host1,inicial_host2):
	##############################
	####CENARIOS COM LIBERAÇÃO####
	##############################
	#liberando 6k cotesia e 7k tetrastichus por hectare
	if loop==29: #100h1_0h2_q
		Pnumerico=52.5
		Q1numerico=30
		lin=50
		col=50 
		replicas=100 #replicas
		TS=[4,4,4] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,3.07]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=8.01
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=1.0341
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.125,0.125,0.05,0.05]#era o antigo com ele migrando para 8 vizinhos#PQNH
		tempo_final=200

		inicial_esp=5
		inicial_gen=5
		#
		solturas=[1]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[5],[5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[0.9,0.9,0.3,0.3]#taxa de qnts parasitoides saem do patch baixa e alta. 
		#a porcentagem vai mudar entre as simulações, e eu preciso seguir a mesma ordem para as duas sp.
		porcentagens_n=[0.0004] #ocupacao inicial do host primario
		porcentagens_h=[0.0004] #ocupacao inicial do host alternativo
		tempo=condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		#wd_novo=f'{wd_original}\\simulacoes'
		#os.makedirs(wd_novo)
		#os.chdir(wd_novo)
		for zzzz in range(replicas): 
			modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico)
	if loop==30: #50h1_50h2_q
		Pnumerico=52.5
		Q1numerico=30
		lin=50
		col=50 
		replicas=100 #replicas
		TS=[4,4,4] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,3.07]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=8.01
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=1.0341
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.125,0.125,0.05,0.05]#era o antigo com ele migrando para 8 vizinhos#PQNH
		tempo_final=200

		inicial_esp=5
		inicial_gen=5
		#
		solturas=[1]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[5],[5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[0.9,0.9,0.3,0.3]#taxa de qnts parasitoides saem do patch baixa e alta. 
		#a porcentagem vai mudar entre as simulações, e eu preciso seguir a mesma ordem para as duas sp.
		porcentagens_n=[0.0004] #ocupacao inicial do host primario
		porcentagens_h=[0.0004] #ocupacao inicial do host alternativo
		tempo=condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		#wd_novo=f'{wd_original}\\simulacoes'
		#os.makedirs(wd_novo)
		#os.chdir(wd_novo)
		for zzzz in range(replicas): 
			modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico)
	if loop==31: #75h1_25h2_q
		Pnumerico=52.5
		Q1numerico=30
		lin=50
		col=50 
		replicas=100 #replicas
		TS=[4,4,4] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,3.07]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=8.01
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=1.0341
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.125,0.125,0.05,0.05]#era o antigo com ele migrando para 8 vizinhos#PQNH
		tempo_final=200

		inicial_esp=5
		inicial_gen=5
		#
		solturas=[1]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[5],[5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[0.9,0.9,0.3,0.3]#taxa de qnts parasitoides saem do patch baixa e alta. 
		#a porcentagem vai mudar entre as simulações, e eu preciso seguir a mesma ordem para as duas sp.
		porcentagens_n=[0.0004] #ocupacao inicial do host primario
		porcentagens_h=[0.0004] #ocupacao inicial do host alternativo
		tempo=condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		#wd_novo=f'{wd_original}\\simulacoes'
		#os.makedirs(wd_novo)
		#os.chdir(wd_novo)
		for zzzz in range(replicas): 
			modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico)
	if loop==32: #25h1_75h2_q
		Pnumerico=52.5
		Q1numerico=30
		lin=50
		col=50 
		replicas=100 #replicas
		TS=[4,4,4] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,3.07]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=8.01
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=1.0341
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.125,0.125,0.05,0.05]#era o antigo com ele migrando para 8 vizinhos#PQNH
		tempo_final=200

		inicial_esp=5
		inicial_gen=5
		#
		solturas=[1]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[5],[5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[0.9,0.9,0.3,0.3]#taxa de qnts parasitoides saem do patch baixa e alta. 
		#a porcentagem vai mudar entre as simulações, e eu preciso seguir a mesma ordem para as duas sp.
		porcentagens_n=[0.0004] #ocupacao inicial do host primario
		porcentagens_h=[0.0004] #ocupacao inicial do host alternativo
		tempo=condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		#wd_novo=f'{wd_original}\\simulacoes'
		#os.makedirs(wd_novo)
		#os.chdir(wd_novo)
		for zzzz in range(replicas): 
			modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico)	
	if loop==33: #100%_100%
		Pnumerico=52.5
		Q1numerico=30
		lin=50
		col=50 
		replicas=100 #replicas
		TS=[4,4,4] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,3.07]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=8.01
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=1.0341
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.125,0.125,0.05,0.05]#era o antigo com ele migrando para 8 vizinhos#PQNH
		tempo_final=200

		inicial_esp=5
		inicial_gen=5
		#
		solturas=[1]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[5],[5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[0.9,0.9,0.3,0.3]#taxa de qnts parasitoides saem do patch baixa e alta. 
		#a porcentagem vai mudar entre as simulações, e eu preciso seguir a mesma ordem para as duas sp.
		porcentagens_n=[0.0004] #ocupacao inicial do host primario
		porcentagens_h=[0.0004] #ocupacao inicial do host alternativo
		tempo=condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		#wd_novo=f'{wd_original}\\simulacoes'
		#os.makedirs(wd_novo)
		#os.chdir(wd_novo)
		for zzzz in range(replicas): 
			modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico)		
	if loop==34: #75h1_25h2_b
		Pnumerico=52.5
		Q1numerico=30
		lin=50
		col=50 
		replicas=100 #replicas
		TS=[4,4,4] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,3.07]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=8.01
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=1.0341
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.125,0.125,0.05,0.05]#era o antigo com ele migrando para 8 vizinhos#PQNH
		tempo_final=200

		inicial_esp=5
		inicial_gen=5
		#
		solturas=[1]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[5],[5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[0.9,0.9,0.3,0.3]#taxa de qnts parasitoides saem do patch baixa e alta. 
		#a porcentagem vai mudar entre as simulações, e eu preciso seguir a mesma ordem para as duas sp.
		porcentagens_n=[0.0004] #ocupacao inicial do host primario
		porcentagens_h=[0.0004] #ocupacao inicial do host alternativo
		tempo=condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		#wd_novo=f'{wd_original}\\simulacoes'
		#os.makedirs(wd_novo)
		#os.chdir(wd_novo)
		for zzzz in range(replicas): 
			modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico)			
	if loop==35: #50h1_50h2_b
		Pnumerico=52.5
		Q1numerico=30
		lin=50
		col=50 
		replicas=100 #replicas
		TS=[4,4,4] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,3.07]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=8.01
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=1.0341
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.125,0.125,0.05,0.05]#era o antigo com ele migrando para 8 vizinhos#PQNH
		tempo_final=200

		inicial_esp=5
		inicial_gen=5
		#
		solturas=[1]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[5],[5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[0.9,0.9,0.3,0.3]#taxa de qnts parasitoides saem do patch baixa e alta. 
		#a porcentagem vai mudar entre as simulações, e eu preciso seguir a mesma ordem para as duas sp.
		porcentagens_n=[0.0004] #ocupacao inicial do host primario
		porcentagens_h=[0.0004] #ocupacao inicial do host alternativo
		tempo=condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		#wd_novo=f'{wd_original}\\simulacoes'
		#os.makedirs(wd_novo)
		#os.chdir(wd_novo)
		for zzzz in range(replicas): 
			modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico)				
	if loop==36: #25h1_75h2_b
		Pnumerico=52.5
		Q1numerico=30
		lin=50
		col=50 
		replicas=100 #replicas
		TS=[4,4,4] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,3.07]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=8.01
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=1.0341
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.125,0.125,0.05,0.05]#era o antigo com ele migrando para 8 vizinhos#PQNH
		tempo_final=200

		inicial_esp=5
		inicial_gen=5
		#
		solturas=[1]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[5],[5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[0.9,0.9,0.3,0.3]#taxa de qnts parasitoides saem do patch baixa e alta. 
		#a porcentagem vai mudar entre as simulações, e eu preciso seguir a mesma ordem para as duas sp.
		porcentagens_n=[0.0004] #ocupacao inicial do host primario
		porcentagens_h=[0.0004] #ocupacao inicial do host alternativo
		tempo=condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		#wd_novo=f'{wd_original}\\simulacoes'
		#os.makedirs(wd_novo)
		#os.chdir(wd_novo)
		for zzzz in range(replicas): 
			modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico)
	if loop==37: #50h1_50h2_q2
		Pnumerico=52.5
		Q1numerico=30
		lin=50
		col=50 
		replicas=100 #replicas
		TS=[4,4,4] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,3.07]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=8.01
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=1.0341
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.125,0.125,0.05,0.05]#era o antigo com ele migrando para 8 vizinhos#PQNH
		tempo_final=200

		inicial_esp=5
		inicial_gen=5
		#
		solturas=[1]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[5],[5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[0.9,0.9,0.3,0.3]#taxa de qnts parasitoides saem do patch baixa e alta. 
		#a porcentagem vai mudar entre as simulações, e eu preciso seguir a mesma ordem para as duas sp.
		porcentagens_n=[0.0004] #ocupacao inicial do host primario
		porcentagens_h=[0.0004] #ocupacao inicial do host alternativo
		tempo=condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		#wd_novo=f'{wd_original}\\simulacoes'
		#os.makedirs(wd_novo)
		#os.chdir(wd_novo)
		for zzzz in range(replicas): 
			modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico)							
	if loop in [38,39,40,41,42,43,44]: #o h1 sempre 100%: 100h1_75h2q,100h1_50h2q,100h1_50h2q2,100h1_25h2q,100h1_25h2b,100h1_50h2b,100h1_75h2b
		Pnumerico=52.5
		Q1numerico=30
		lin=50
		col=50 
		replicas=100 #replicas
		TS=[4,4,4] #esse valor é o que a outra aluna da professora achou
		TH=[1.46,3.07,3.07]#esse valor é o que a outra aluna da professora achou
		a1=0.28 #esse valor é o que a outra aluna da professora achou
		a2=8.01 #esse valor é o que a outra aluna da professora achou
		a3=8.01
		lambda1=1.0341#esse é o valor tirado da tese da tabela de vida desse bixo
		lambda2=1.0341
		h0=[100,100]
		H0=[300,300]
		f0=[100,100]
		taxa_disp=[0.125,0.125,0.05,0.05]#era o antigo com ele migrando para 8 vizinhos#PQNH
		tempo_final=200

		inicial_esp=5
		inicial_gen=5
		#
		solturas=[1]
		wd_original=f'{os.getcwd()}'
		para_lista_inicial=[[5],[5]]#to-do: acomodar esse novo jeito nos outros loops e na função
		
		
		lista_raios=viz[1] #salvando a lista de raios
		viz=viz[0] #salvando a lista de vizinhos
		indv_migrante=[0.9,0.9,0.3,0.3]#taxa de qnts parasitoides saem do patch baixa e alta. 
		#a porcentagem vai mudar entre as simulações, e eu preciso seguir a mesma ordem para as duas sp.
		porcentagens_n=[0.0004] #ocupacao inicial do host primario
		porcentagens_h=[0.0004] #ocupacao inicial do host alternativo
		tempo=condicao_inicial_paisagem_infesstada_100_parasitoide_meio(lin,col,inicial_host1,inicial_host2,inicial_esp,inicial_gen,loop)
		patch_iniciais=tempo[0]
		pop_iniciais=tempo[1]
		#wd_novo=f'{wd_original}\\simulacoes'
		#os.makedirs(wd_novo)
		#os.chdir(wd_novo)
		for zzzz in range(replicas): 
			modelo_paisagem_migracao(a1,a2,a3,lambda1,lambda2,H0,h0,f0,taxa_disp,indv_migrante,tempo_final,lin,col,pop_iniciais,patch_iniciais,viz,lista_raios,TS,TH,zzzz,para_lista_inicial,loop,lib_extras,Pnumerico,Q1numerico)								
	
	
	return






############################
############################
############################


lib_extras=False
loop=29
lin=50
col=50
distancia_maxima=3
lista_distancia=checando_existencia_lista_distancias(lin,col,distancia_maxima)
tempoo=checando_existencia_lista_viz_e_lista_de_patches(lin,col,distancia_maxima)
viz=tempoo[0]
lista_patches=tempoo[1]
lista_borda=tempoo[2]
contagem_especialista50=False
contagem_especialista100=False
contagem_especialista200=False

if loop==29:
	wd_1=f'{os.getcwd()}'
	os.chdir(wd_1)
	for loop_numero in [29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44]:
		if loop_numero==29: #100 e 0
			os.chdir(wd_1)
			inicial_host1=100
			inicial_host2=0
			pasta=f'100h1_0h2_q'
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==30:
			inicial_host1=100
			inicial_host2=100
			pasta=f'50h1_50h2_q'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==31:
			inicial_host1=100 #mas 75% da paisagem
			inicial_host2=100 #mas 25% da paisagem
			pasta=f'75h1_25h2_q'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==32:
			inicial_host1=100 #mas 25% da paisagem
			inicial_host2=100 #mas 75% da paisagem
			pasta=f'25h1_75h2_q'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==33:
			inicial_host1=100 #mas 100% da paisagem
			inicial_host2=100 #mas 100% da paisagem
			pasta=f'100%_100%'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==34:
			inicial_host1=100 #mas 73.92% da paisagem
			inicial_host2=100 #mas 26.08% da paisagem só q encostado em todas as bordas
			pasta=f'75h1_25h2_b'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==35:
			inicial_host1=100 #50%
			inicial_host2=100 #50%
			pasta=f'50h1_50h2_b'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==36:
			inicial_host1=100 #25%
			inicial_host2=100 #75%
			pasta=f'25h1_75h2_b'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==37:
			inicial_host1=100 #50%
			inicial_host2=100 #50%
			pasta=f'50h1_50h2_q2'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==38:
			inicial_host1=100
			inicial_host2=100
			pasta=f'100h1_75h2_q'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==39:
			inicial_host1=100 #mas 75% da paisagem
			inicial_host2=100 #mas 25% da paisagem
			pasta=f'100h1_50h2_q'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==40:
			inicial_host1=100 #mas 25% da paisagem
			inicial_host2=100 #mas 75% da paisagem
			pasta=f'100h1_50h2q_2'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==41:
			inicial_host1=100 #mas 75% da paisagem
			inicial_host2=100 #mas 25% da paisagem só q numa faixa do canto
			pasta=f'100h1_25h2_q'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==42:
			inicial_host1=100 #mas 73.92% da paisagem
			inicial_host2=100 #mas 26.08% da paisagem só q encostado em todas as bordas
			pasta=f'100h1_25h2_b'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==43:
			inicial_host1=100 #50%
			inicial_host2=100 #50%
			pasta=f'100h1_50h2_b'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		if loop_numero==44:
			inicial_host1=100 #25%
			inicial_host2=100 #75%
			pasta=f'100h1_75h2_b'
			os.chdir(wd_1)
			os.makedirs(pasta)
			os.chdir(pasta)
		loops_cenarios(loop_numero,viz,lista_patches,lista_borda,inicial_host1,inicial_host2)		



####
#oque falta mudar:
#modelo_paisagem_migraçao: ok
#criando_lista_patches_cada_hectare: ok
#criando_lista_patches_cada_borda: ok
#criando_lista_distancia: ok
#checando_existencia_lista_distancias: ok
#checando_existencia_lista_viz_e_lista_de_patches: ok
#criando_lista_viz: ok
#criando_popmigrante_densidade_patch: Ok
#sorteio_vizinhos: ok
#dinamica_tempo0: Ok #to-do: parece estar ok, mas dependendo das condições iniciais, eu tenha q mudar, mas acho q mesmo assim ja esta ok
#dinamica_iterativa: OK
#conferir_recurso_parasitoides: Ok
#salva_arquivos: ok
#criador_vetor: ok
#sorteio_pop_inicial_host: OK #to-do: falta arrumar os loops de cenarios para corresponder aos cenarios q vamos usar
#alocar_parasitoide: #to-do: falta arrumar os loops de cenarios para corresponder aos cenarios q vamo
#MEIO_alocando_parasitoide: #to-do: esse tem que mudar quando for definir as simulações
#BORDA_alocando_parasitoide: #to-do: esse tem que mudar quando for definir as simulações
#checando_necessidade_liberacao: --
#nova_liberacao: --
#loops_cenarios: #to-do: falta arrumar os loops de cenarios para corresponder aos cenarios q vamos 