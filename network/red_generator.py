import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
import igraph

# class Red():
#     '''
#         Attributes : -> Q = total number of cliques (int)
#                     -> n_max = maximum size a cliue can have (int)
#                     -> gamma = probability of a node having an edge to connecting to other clique (float)
#                     -> fn = probability distribution for a clique being of size n (function)
#                     -> clique_list = clique list (list[nx.Graph])
#                     -> red_total = full network (nx.Graph)
#     '''

#     def __init__(self, Q, Gamma , fn, n_max):
                
#         # Asigno valores y creo la lista de cliques:
#         self.Q = Q  # cantidad total de cliques
#         self.gamma = Gamma # probabilidad de que un nodo tenga enlace a otro clique
#         self.n_max = n_max # maximo tamaño de clique
#         self.fn = fn  # distribucion de probabilidad para el tamaño de clique
#         self.clique_list = [] # lista con los cliques de la red 
        

#         # Creo una lista de nombres de la forma "clique 'nro_clique' nodo 'nro_nodo'" :
#         name_list = []

#         # Genero una seed aleatoria:
#         np.random.seed()

#         # Creo un array de los posibles valores de tamaño de clique :
#         n_vals = np.arange(1, n_max + 1)
#         f_n_vals = fn(n_vals)
#         f_n_vals = f_n_vals / np.sum(f_n_vals) # normalizo

#         for i in range(Q):
#             name_list.append("clique " + str(i) + " nodo ")
            
#             n = np.random.choice(n_vals, p = f_n_vals)
            
#             aux_clique = igraph.Graph.Full(n = n, directed=False, loops = False) # genero el clique
#             self.clique_list.append(aux_clique)

#     #     Pongo que los nodos tengan como atributo el numero de clique al que pertenecen:
#         cuenta = 0
#         for clique in self.clique_list:
#             for  node in clique.vs.indices:
#                 clique.vs[node]["nro_clique"] = cuenta
#                 clique.vs[node]["ext_connection"] = False

#             cuenta += 1        
        
#         self.red_total_cliques_desconectados = igraph.disjoint_union(self.clique_list) # red_total es un grafo de la red completa
        
         
#         self.red_total_cliques_conectados = self.red_total_cliques_desconectados.copy()

#     #     # # Asigno los enlaces entre cliques :
#         aux_node_list = []
#         for nodo in self.red_total_cliques_conectados.vs.indices:
#             nro_random = np.random.choice([0,1], p = [1 - self.gamma, self.gamma])        
#             if nro_random:
#                 aux_node_list.append(nodo)
            
#         np.random.shuffle(aux_node_list)

#         while len(aux_node_list) > 1:
            
            
#             nodo0 = aux_node_list[0]
            
#             l = len(aux_node_list)

#             clique_0 = self.red_total_cliques_conectados.vs[nodo0]["nro_clique"]
#             cuenta = 0
#             clique_cuenta = self.red_total_cliques_conectados.vs[nodo0]["nro_clique"]
#             while clique_cuenta == clique_0:
#                 if cuenta == l:
#                     break
#                 clique_cuenta = self.red_total_cliques_conectados.vs[aux_node_list[cuenta]]["nro_clique"]
#                 cuenta += 1
                
#             if cuenta >= l:
#                 break
            
#             nodo1 = aux_node_list[1]

#             while self.red_total_cliques_conectados.are_connected(nodo0, nodo1):
#                 np.random.shuffle(aux_node_list)
#                 nodo0 = aux_node_list[0]
#                 nodo1 = aux_node_list[1]         

#             self.red_total_cliques_conectados.add_edges([(nodo0, nodo1)])
#             aux_node_list.remove(nodo0)
#             aux_node_list.remove(nodo1)
#             self.red_total_cliques_conectados.vs[nodo0]["ext_connection"] = True
#             self.red_total_cliques_conectados.vs[nodo1]["ext_connection"] = True

#     def guardar_archivo_graphml(self, file_name):
#         red.red_total_cliques_conectados.write_graphml(file_name + ".graphml")

#     def grado(self, nodo):
#         return self.red_total_cliques_conectados.vs[nodo].degree()

#     def clustering_local(self, nodo):
#         return self.red_total_cliques_conectados.transitivity_local_undirected(self.red_total_cliques_conectados.vs[nodo], mode = "zero")

#     def clustering_medio(self):
#         return self.red_total_cliques_conectados.transitivity_avglocal_undirected(mode = 'zero')
    
#     def clustering_global(self):
#         return self.red_total_cliques_conectados.transitivity_undirected(mode = "zero")

#     def nodos(self):
#         return self.red_total_cliques_conectados.vs.indices

#     def cantidad_de_nodos(self):
#         return len(self.nodos())

#     def enlaces(self):

#         ret_list = []
#         for edge in self.red_total_cliques_conectados.es:
#             source_vertex_id = edge.source
#             target_vertex_id = edge.target
#             source_vertex = self.red_total_cliques_conectados.vs.indices[source_vertex_id]
#             target_vertex = self.red_total_cliques_conectados.vs.indices[target_vertex_id]
#             ret_list.append((source_vertex, target_vertex))

#         return ret_list

#     def cantidad_de_enlaces(self):
#         return len(self.enlaces())

#     def gamma_efectivo(self):
#         num_nodos_con_enlace_externo = sum(2 for enlace in self.enlaces() if self.red_total_cliques_conectados.vs[enlace[0]]["nro_clique"] != self.red_total_cliques_conectados.vs[enlace[1]]["nro_clique"])
#         return num_nodos_con_enlace_externo / len(self.nodos())

#     def clustering_medio_componente_mas_grande(self):
#         CG = self.red_total_cliques_conectados.connected_components().giant()
#         return CG.transitivity_avglocal_undirected(mode = 'zero')

#     def clustering_global_componente_mas_grande(self):
#         CG = self.red_total_cliques_conectados.connected_components().giant()
#         return CG.transitivity_undirected(mode = "zero")
    
#     def mayor_componente(self):
#         CG = self.red_total_cliques_conectados.connected_components().giant()
#         return CG

#     def cantidad_de_componentes(self):
#         return len(list(self.red_total_cliques_conectados.connected_components()))

#     def cantidad_de_cliques_en_mayor_componente(self):
#         CG = self.red_total_cliques_conectados.connected_components().giant()

#         clique_number_list = [CG.vs[nodo]["nro_clique"] for nodo in CG.vs.indices]
#         clique_number_set = set(clique_number_list)
            
#         return len(clique_number_set)

#     def cuenta_cliques_sizes_en_componente_mas_grande(self):
        
#         CG = self.red_total_cliques_conectados.connected_components().giant()

#         clique_number_list = [CG.vs[nodo]["nro_clique"] for nodo in CG.vs.indices]
#         clique_number_list = list(set(clique_number_list))  

#         sizes = np.zeros(self.n_max)
#         for i in range(len(clique_number_list)):
#             nro_clique = clique_number_list[i]
#             size = len(self.clique_list[nro_clique].vs.indices)
#             sizes[size-1] += 1

#         return sizes


#     # def average_path_componente_gigante(self):
#     #     Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
#     #     componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])

#     #     return nx.average_shortest_path_length(componente_gigante)

#     # # def closeness_centrality(self, nodo):
#     # #     return nx.closeness_centrality(self.red_total_cliques_conectados, u = nodo)



#     # def cuenta_cliques_sizes(self):
#     #     sizes = np.zeros(self.n_max)
#     #     for i in range(len(self.clique_list)):
#     #         size = len(self.clique_list[i])
#     #         sizes[size-1] += 1

#     #     return sizes        
    
#     # def homofilia_grado(self):
#     #     return nx.degree_assortativity_coefficient(self.red_total_cliques_conectados)

#     # def homofilia_grado_componente_gigante(self):
#     #     Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
#     #     componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])
#     #     return nx.degree_assortativity_coefficient(componente_gigante)

#     # def fraccion_nodos_elegibles_en_CG(self):
#     #     Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
#     #     componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])
#     #     enlace_externo = nx.get_node_attributes(componente_gigante, "ext_connection")
#     #     cuenta = 0
#     #     for nodo in componente_gigante.nodes():
#     #         if not enlace_externo[nodo]:
#     #             cuenta += 1

#     #     return cuenta / len(self.nodes())




class Red():
    '''
        Attributes : -> Q = total number of cliques (int)
                    -> n_max = maximum size a cliue can have (int)
                    -> gamma = probability of a node having an edge to connecting to other clique (float)
                    -> fn = probability distribution for a clique being of size n (function)
                    -> clique_list = clique list (list[nx.Graph])
                    -> red_total = full network (nx.Graph)
    '''

    def __init__(self, Q, Gamma , fn, n_max):
                
        # Asigno valores y creo la lista de cliques:
        self.Q = Q  # cantidad total de cliques
        self.gamma = Gamma # probabilidad de que un nodo tenga enlace a otro clique
        self.n_max = n_max # maximo tamaño de clique
        self.fn = fn  # distribucion de probabilidad para el tamaño de clique
        self.clique_list = [] # lista con los cliques de la red 
        

        # Genero una seed aleatoria:
        np.random.seed()

        # Creo un array de los posibles valores de tamaño de clique :
        n_vals = np.arange(1, n_max + 1)
        f_n_vals = fn(n_vals)
        f_n_vals = f_n_vals / np.sum(f_n_vals) # normalizo

        for i in range(Q):            
            n = np.random.choice(n_vals, p = f_n_vals)
            
            aux_clique = nx.complete_graph(n) # genero el clique
            self.clique_list.append(aux_clique)

        cuenta = 0
        for clique in self.clique_list:
            for  node in clique.nodes:
                clique.nodes[node]["nro_clique"] = cuenta
                clique.nodes[node]["ext_connection"] = False

            cuenta += 1  

        self.red_total_cliques_desconectados = nx.disjoint_union_all(self.clique_list) #, rename = name_list) # red_total es un grafo de la red completa
        nx.set_edge_attributes(self.red_total_cliques_desconectados, False, name = "largo_alcance")

        self.red_total_cliques_conectados = self.red_total_cliques_desconectados.copy()

        # # Asigno los enlaces entre cliques :
        
        aux_node_list = []
        for nodo in self.red_total_cliques_conectados.nodes():
            nro_random = np.random.choice([0,1], p = [1 - self.gamma, self.gamma])        
            if nro_random:
                aux_node_list.append(nodo)
            
        np.random.shuffle(aux_node_list)

        while len(aux_node_list) > 1:
            
            
            nodo0 = aux_node_list[0]
            
            l = len(aux_node_list)

            clique_0 = self.red_total_cliques_conectados.nodes[nodo0]["nro_clique"]
            cuenta = 0
            clique_cuenta = self.red_total_cliques_conectados.nodes[nodo0]["nro_clique"]
            while clique_cuenta == clique_0:
                if cuenta == l:
                    break
                clique_cuenta = self.red_total_cliques_conectados.nodes[aux_node_list[cuenta]]["nro_clique"]
                cuenta += 1
                
            if cuenta >= l:
                break
            
            nodo1 = aux_node_list[1]

            while self.red_total_cliques_conectados.has_edge(nodo0, nodo1):
                np.random.shuffle(aux_node_list)
                nodo0 = aux_node_list[0]
                nodo1 = aux_node_list[1]         

            self.red_total_cliques_conectados.add_edge(nodo0, nodo1, largo_alcance = True)

            aux_node_list.remove(nodo0)
            aux_node_list.remove(nodo1)
            self.red_total_cliques_conectados.nodes[nodo0]["ext_connection"] = True
            self.red_total_cliques_conectados.nodes[nodo1]["ext_connection"] = True

        
         


    def guardar_archivo_graphml(self, file_name):
        grafo = igraph.Graph.from_networkx(self.red_total_cliques_conectados)
        grafo.write_graphml(file_name + ".graphml")

    # def grado(self, nodo):
    #     return self.red_total_cliques_conectados.vs[nodo].degree()

    # def clustering_local(self, nodo):
    #     return self.red_total_cliques_conectados.transitivity_local_undirected(self.red_total_cliques_conectados.vs[nodo], mode = "zero")

    # def clustering_medio(self):
    #     return self.red_total_cliques_conectados.transitivity_avglocal_undirected(mode = 'zero')
    
    # def clustering_global(self):
    #     return self.red_total_cliques_conectados.transitivity_undirected(mode = "zero")

    # def nodos(self):
    #     return self.red_total_cliques_conectados.vs.indices

    # def cantidad_de_nodos(self):
    #     return len(self.nodos())

    # def enlaces(self):

    #     ret_list = []
    #     for edge in self.red_total_cliques_conectados.es:
    #         source_vertex_id = edge.source
    #         target_vertex_id = edge.target
    #         source_vertex = self.red_total_cliques_conectados.vs.indices[source_vertex_id]
    #         target_vertex = self.red_total_cliques_conectados.vs.indices[target_vertex_id]
    #         ret_list.append((source_vertex, target_vertex))

    #     return ret_list

    # def cantidad_de_enlaces(self):
    #     return len(self.enlaces())

    # def gamma_efectivo(self):
    #     num_nodos_con_enlace_externo = sum(2 for enlace in self.enlaces() if self.red_total_cliques_conectados.vs[enlace[0]]["nro_clique"] != self.red_total_cliques_conectados.vs[enlace[1]]["nro_clique"])
    #     return num_nodos_con_enlace_externo / len(self.nodos())

    # def clustering_medio_componente_mas_grande(self):
    #     CG = self.red_total_cliques_conectados.connected_components().giant()
    #     return CG.transitivity_avglocal_undirected(mode = 'zero')

    # def clustering_global_componente_mas_grande(self):
    #     CG = self.red_total_cliques_conectados.connected_components().giant()
    #     return CG.transitivity_undirected(mode = "zero")
    
    # def mayor_componente(self):
    #     CG = self.red_total_cliques_conectados.connected_components().giant()
    #     return CG

    # def cantidad_de_componentes(self):
    #     return len(list(self.red_total_cliques_conectados.connected_components()))

    # def cantidad_de_cliques_en_mayor_componente(self):
    #     CG = self.red_total_cliques_conectados.connected_components().giant()

    #     clique_number_list = [CG.vs[nodo]["nro_clique"] for nodo in CG.vs.indices]
    #     clique_number_set = set(clique_number_list)
            
    #     return len(clique_number_set)

    # def cuenta_cliques_sizes_en_componente_mas_grande(self):
        
    #     CG = self.red_total_cliques_conectados.connected_components().giant()

    #     clique_number_list = [CG.vs[nodo]["nro_clique"] for nodo in CG.vs.indices]
    #     clique_number_list = list(set(clique_number_list))  

    #     sizes = np.zeros(self.n_max)
    #     for i in range(len(clique_number_list)):
    #         nro_clique = clique_number_list[i]
    #         size = len(self.clique_list[nro_clique].vs.indices)
    #         sizes[size-1] += 1

    #     return sizes

######################################## ALGORITMOS CON NETWORKX:

    def plot_conectado(self):
        nx.draw(self.red_total_cliques_conectados)
    
    def plot_desconectado(self):
        nx.draw(self.red_total_cliques_desconectados)
    
    def plot_componente_mas_grande(self):
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])
        nx.draw(componente_gigante)

    def degree(self, node):
        return self.red_total_cliques_conectados.degree[node]

    def distribucion_de_grado_sin_normalizacion(self):
    
        cant_elementos_en_grado_k = np.zeros(self.n_max + 1)
        for nodo in self.red_total_cliques_conectados.nodes():
            grado = self.red_total_cliques_conectados.degree[nodo]
            cant_elementos_en_grado_k[grado] += 1

        return cant_elementos_en_grado_k

    def distribucion_de_grado_sin_normalizacion_componente_gigante(self):
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])
        
        cant_elementos_en_grado_k = np.zeros(self.n_max + 1)
        for nodo in componente_gigante.nodes():
            grado = self.red_total_cliques_conectados.degree[nodo]
            cant_elementos_en_grado_k[grado] += 1

        return cant_elementos_en_grado_k


    def clustering_local(self, nodo):
        return nx.clustering(self.red_total_cliques_conectados, nodo)

    def clustering_medio(self):
        return nx.average_clustering(self.red_total_cliques_conectados, count_zeros = True)

    def clustering_global(self):
        return nx.transitivity(self.red_total_cliques_conectados)
    
    def clustering_global_igraph(self):
        return self.red_total_cliques_conectados_igraph.transitivity_undirected()

    def nodes(self):
        return self.red_total_cliques_conectados.nodes()

    def cantidad_de_nodos(self):
        return len(self.red_total_cliques_conectados.nodes())

    def cantidad_de_nodos_componente_mas_grande(self):
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        largest_cc = self.red_total_cliques_conectados.subgraph(Gcc[0])
        
        return len(largest_cc.nodes())

    def edges(self):
        return self.red_total_cliques_conectados.edges()

    def cantidad_de_enlaces(self):
        return len(self.red_total_cliques_conectados.edges())

    def gamma_efectivo(self):
        num_nodos_con_enlace_externo = sum(2 for edge in self.edges() if self.red_total_cliques_conectados.edges[edge]["largo_alcance"])
        return num_nodos_con_enlace_externo / len(self.nodes())

    def clustering_medio_componente_mas_grande(self):
        
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        largest_cc = self.red_total_cliques_conectados.subgraph(Gcc[0])

        clust_return = nx.average_clustering(largest_cc, count_zeros = True)
        return clust_return

    def clustering_global_componente_mas_grande(self):
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        largest_cc = self.red_total_cliques_conectados.subgraph(Gcc[0])
        clust_return = nx.transitivity(largest_cc)
        return clust_return
    
    def mayor_componente(self):
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        G0 = self.red_total_cliques_conectados.subgraph(Gcc[0])
        return G0

    def cantidad_de_componentes(self):
        return len(list(nx.connected_components(self.red_total_cliques_conectados)))

    def component_sizes(self): 
        Gcc = nx.connected_components(self.red_total_cliques_conectados)
        return [len(g) for g in Gcc]

#     def average_path_componente_gigante(self):
#         Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
#         componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])

#         return nx.average_shortest_path_length(componente_gigante)


    def cantidad_de_cliques_en_mayor_componente(self):
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        G0 = self.red_total_cliques_conectados.subgraph(Gcc[0])


        clique_number_list = [G0.nodes[nodo]["nro_clique"] for nodo in G0.nodes()]
        clique_number_set = set(clique_number_list)
            
        return len(clique_number_set)

    def cuenta_cliques_sizes_en_componente_mas_grande(self):
        
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])

        clique_number_list = [componente_gigante.nodes[nodo]["nro_clique"] for nodo in componente_gigante.nodes()]
        clique_number_list = list(set(clique_number_list))  

        sizes = np.zeros(self.n_max)
        for i in range(len(clique_number_list)):
            nro_clique = clique_number_list[i]
            size = len(self.clique_list[nro_clique])
            sizes[size-1] += 1

        return sizes

    def cuenta_cliques_sizes(self):
        sizes = np.zeros(self.n_max)
        for i in range(len(self.clique_list)):
            size = len(self.clique_list[i])
            sizes[size-1] += 1

        return sizes        
    
    def homofilia_grado(self):
        suma1 = 0
        suma2 = 0
        suma3 = 0
        
        for enlace in self.red_total_cliques_conectados.edges():
            k1 = self.degree(enlace[0])
            k2 = self.degree(enlace[1])

            suma1 += k1 * k2
            suma2 += k1 + k2
            suma3 += k1**2 + k2 **2

        suma1 /= len(self.red_total_cliques_conectados.edges())
        suma2 /= len(self.red_total_cliques_conectados.edges()) 
        suma3 /= len(self.red_total_cliques_conectados.edges())       

        return suma1, suma2, suma3

    def homofilia_grado_componente_gigante(self):
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])

        suma1 = 0
        suma2 = 0
        suma3 = 0
        
        for enlace in componente_gigante.edges():
            k1 = self.degree(enlace[0])
            k2 = self.degree(enlace[1])

            suma1 += k1 * k2
            suma2 += k1 + k2
            suma3 += k1**2 + k2 **2

        suma1 /= len(componente_gigante.edges())
        suma2 /= len(componente_gigante.edges()) 
        suma3 /= len(componente_gigante.edges())       

        return suma1, suma2, suma3


    def fraccion_nodos_con_enlace_CG(self):
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])
        enlace_externo = nx.get_node_attributes(componente_gigante, "ext_connection")
        cuenta = 0
        for nodo in componente_gigante.nodes():
            if enlace_externo[nodo]:
                cuenta += 1

        return cuenta / len(componente_gigante.nodes())    

    def media_menor_distancia(self):
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])

        CG_igraph = igraph.Graph.from_networkx(componente_gigante)

        return CG_igraph.average_path_length()

    def media_menor_distancia_red_total(self):
        red_igraph = igraph.Graph.from_networkx(self.red_total_cliques_conectados)

        return red_igraph.average_path_length(directed = False, unconn = True)


    def max_menor_distancia(self):
        # Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        # componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])
        # return nx.diameter(componente_gigante)
        Gcc = sorted(nx.connected_components(self.red_total_cliques_conectados), key=len, reverse=True)
        componente_gigante = self.red_total_cliques_conectados.subgraph(Gcc[0])

        CG_igraph = igraph.Graph.from_networkx(componente_gigante)

        return CG_igraph.diameter()



if __name__ == "__main__":


    def fn(n):
        return np.exp(-0.05 * n)

    from my_package.my_functions import distributions as mis_distr



    fun_kro = mis_distr.fn_constante(5)

    red = Red(10, 0.5, fn, 10)

    red.guardar_archivo_graphml("C:/Users/Lucas/Documents/IB/MAESTRIA/Figuras trabajo especial/red_ejemplo_numerico_nuevo")
    print(red.gamma_efectivo())

    print("---main end---")