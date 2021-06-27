from cmath import *
import numpy as numpy
import matplotlib.animation as animation
import matplotlib.patheffects as path_effects
import sys
sys.path.append('/Users/dasilva/Google Drive/Thèse/Python/')
from mouvements import *
from math import factorial
import networkx as nx


def f(d):
    if d==1 or d==-1:
        return 0
    else:
        return 1


class Gmod:
    def __init__(self,lpoint,ldirection,lmouv,orienté=False,multi=False):
        self.lpoint=lpoint
        self.ldirection=ldirection
        self.lmouv=lmouv
        if  orienté==False:
            if multi==False:
                Gtemp=nx.Graph()
            else:
                Gtemp=nx.MultiGraph()
        else:
            if multi==False:
                Gtemp=nx.DiGraph()
            else:
                Gtemp=nx.MultiDiGraph()
        N=len(lpoint)
        num=1
        for k in range(N-1):
            Gtemp.add_edge((lpoint[k],f(ldirection[k])),\
                (lpoint[k+1],f(ldirection[k+1])),\
                    num=num,mouv=lmouv[k].__name__)
            num+=1
        self.graphe=Gtemp
        pos=dict() # on attribue des positions aux sommets
        V=list(Gtemp.nodes) #list of vertex
        couleur=[] # liste des couleurs
        for v in V: 
            eps=(-1)**(v[0].imag)
            if v[1]==0:
                pos[v]=(v[0].real+eps*0.1,v[0].imag+0.2)
                couleur+=['r']
            else:
                pos[v]=(v[0].real-eps*0.1,v[0].imag-0.2)
                couleur+=['b']
        self.pos=pos
        self.node_color=couleur

def sonde(L):
    N=len(L)
    for k in range(N):
        for l in range(k+1,N):
            if L[k]==L[l]:
                return (k,l)
    return False


def est_un_cycle(L):
    if L[0]==L[-1] and sonde(L[1:-1])==False:
        return True
    else:
        return False

def decompose_cycle(G): #G est un objet Gmod
    L=[(G.lpoint[k],G.ldirection[k]) for k in range(len(G.lpoint))]
    if est_un_cycle(L):
        return [G]
    else:
        (i,j)=sonde(L[:-1]) # on cherche deux sommets egaux entre 1 et N-2
        #i+=1 # décalage indice
        #j+=1
        G1=Gmod(G.lpoint[i:j+1],G.ldirection[i:j+1],G.lmouv[i:j])
        lpoint2=G.lpoint[:i+1]+G.lpoint[j+1:]
        ldirection2=G.ldirection[:i+1]+G.ldirection[j+1:]
        lmouv2=G.lmouv[:i]+G.lmouv[j:]
        G2=Gmod(lpoint2,ldirection2,lmouv2)
        return decompose_cycle(G1)+decompose_cycle(G2)

def decompose_liste_cycle(L,lmouv): #L = [(lpoint,ldirection)]
    if est_un_cycle(L):
        return [L,lmouv]
    else:
        (i,j)=sonde(L[:-1]) # on cherche deux sommets egaux entre 1 et N-2
        #i+=1 # décalage indice
        #j+=1
        L1=L[:i+1]+L[j+1:]
        lmouv1=lmouv[:i]+lmouv[j:]
        L2=L[i:j+1]
        lmouv2=lmouv[i:j]
        return decompose_liste_cycle(L1,lmouv1)+decompose_liste_cycle(L2,lmouv2)

def coeff(G,u,v):# G un Gmod orienté, u et v deux sommets
        if u==v:
            return G.graphe.out_degree(v)
        else:
            return -1*G.graphe.number_of_edges(u,v)


# retourne le nombre de chemin eulérien par formule de  BEST pour G orienté
def nombre_circuit(G): #G un Gmod
    P=1
    L=list(nx.nodes(G.graphe))
    ######Produit dans la formule de BEST
    for v in L :
        P=P*factorial(G.graphe.out_degree(v)-1) 
        # Compiler avec ça la matrix Laplacienne, calculer n'importe quel mineur,
        # Multiplier ce nombre par P : on obtient le nombre de circuit.
    #M=numpy.matrix([[coeff(u,v) for u in L] for v in L]) Matrice Laplacienne
    Mineur=numpy.matrix([[coeff(G,u,v) for u in L[:-1]] for v in L[:-1]]) #Mineur
    #On enlève la dernière ligne/colonne
    T=numpy.linalg.det(Mineur) #Nombre d'arborescence
    return(int(T*P/(len(list(G.graphe.nodes)))))
      
                     
def encadre_chemin(G): #G Gmod non orienté
    V=list(nx.nodes(G.graphe))
    p=len(V)
    q=len(list(nx.edges(G.graphe)))
    m,M=2**(p-q),2**(p-q)
    ######Produit dans la formule de BEST
    for v in V:
        m=m*factorial((G.graphe.degree(v))/2-1)
        M=M*factorial(G.graphe.degree(v)-1)/factorial((G.graphe.degree(v))/2-1)
    return("Min = "+str(m)+" et Max = "+str(M))