import matplotlib.pyplot as plt
from cmath import *
import numpy as numpy
import sys
sys.path.append('/Users/dasilva/Google Drive/Thèse/Python/')
import matplotlib.patheffects as path_effects
import config #Pour les constantes globales

######Coefficient binomial rapide ########@
def binomial(n, k):
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0


#### Courbe de bézier pour une liste L de complexes, fonction de t #######
def bezier(L,t):
    N=len(L)-1
    return sum([binomial(N,k)*(t**k)*((1-t)**(N-k))*L[k] for k in range(N+1)])


# Segment : seg
def seg(a,b,t):
   return (1-t)*a+t*b

# Virage gauche : vg
def vg(a,b,t):
    z0=a+(b-a)*exp(1j*pi/4)/sqrt(2)
    return z0+exp(1j*pi*t/2)*(a-z0)

def v_2_g(a,b,t):
    z0=a+(b-a)*exp(1j*pi/4)/sqrt(2)
    return z0+exp(1j*pi*t/2)*(a-z0)

# Virage droit : vd
def vd(a,b,t):
    z0=a+(b-a)*exp(-1j*pi/4)/sqrt(2)
    return z0+exp(-1j*pi*t/2)*(a-z0)

def v_2_d(a,b,t):
    z0=a+(b-a)*exp(-1j*pi/4)/sqrt(2)
    return z0+exp(-1j*pi*t/2)*(a-z0)
    
#virage en v à droite
def vvd(a,b,t):
    if t<=0.5:return seg(a,(1j*a-b)/(1j-1),2*t)
    else:return seg((1j*a-b)/(1j-1),b,2*t-1)

def vvg(b,a,t):
    return vvd(a,b,1-t)




 #############################################################@   
#grand virage en v à droite : idem vvd mais se déplace de 2 sur la grille
def Gvvd(a,b,t):
    if t<=0.5:return seg(a,(1j*a-b)/(1j-1),2*t)
    else:return seg((1j*a-b)/(1j-1),b,2*t-1)

def Gvvg(b,a,t):
    return Gvvd(a,b,1-t)
################################################################




# Grand virage gauche : Gvg
def Gvg(a,b,t):
    k=config.taille_Gv
    return ((1 - t)**4)*a + 4*((1 - t)**3)*t*(a + (k)*(b-a)*exp(-1j*3*pi/4))+\
           6*((1 - t)**2)*(t**2)*((a+b)/2-1j*(k)*(b-a)*1.2) + \
           4*(1 - t)*(t**3)*(b +(k)*(a-b)*exp(1j*3*pi/4)) + (t**4)*b

#Grand virage droit : GVd
def Gvd(a,b,t):
    k=config.taille_Gv
    return ((1 - t)**4)*a + 4*((1 - t)**3)*t*(a + (k)*(b-a)*exp(1j*3*pi/4))+\
           6*((1 - t)**2)*(t**2)*((a+b)/2+1j*(k)*(b-a)*1.2) + \
           4*(1 - t)*(t**3)*(b +(k)*(a-b)*exp(-1j*3*pi/4)) + (t**4)*b


# Cercle direct entre a et b
def cercleg(a,b,t):
    z0=a+(b-a)*exp(1j*pi/4)/sqrt(2)
    return z0+exp(2j*pi*t)*(a-z0)

#Cercle indirect entre a et b

def cercled(a,b,t):
    z0=a+(b-a)*exp(-1j*pi/4)/sqrt(2)
    return z0+exp(-2j*pi*t)*(a-z0)



#Petite boucle gauche, point de départ=z, dans la direction du vecteur u (qui donne
# la direction de la boucle et son axe de symétrie. C'est une Bézier-5 points.
# taille_boucle controle la largeur de la boucle


def pbg(z,u,t):
    k=config.taille_boucle
    return ((1 - t)**4)*z + 4*((1 - t)**3)*t*(z+(k*u)/(1+1j))+\
           6*((1 - t)**2)*(t**2)*(z+k*u) + \
           4*(1 - t)*(t**3)*(z+(k*u)/(1-1j)) + (t**4)*z

def pbg_concave(z,u,t):
    k=config.taille_boucle
    return ((1 - t)**4)*z + 4*((1 - t)**3)*t*(z+(0.5*u)/(1+1j))+\
           6*((1 - t)**2)*(t**2)*(z-k*u) + \
           4*(1 - t)*(t**3)*(z+(0.5*u)/(1-1j)) + (t**4)*z

def pbd(z,u,t):
    k=config.taille_boucle
    return ((1 - t)**4)*z + 4*((1 - t)**3)*t*(z+(k*u)/(1-1j))+\
           6*((1 - t)**2)*(t**2)*(z+k*u) + \
           4*(1 - t)*(t**3)*(z+(k*u)/(1+1j)) + (t**4)*z

def pbd_concave(z,u,t):
    return pbg_concave(z,u,1-t)

def demcerclemoins(a,b,t):
    return 1/2*(a+b+(a-b)*exp(-1j*pi*t))

#demi-cercle direct entre a et b
def demcercleplus(a,b,t):
    return 1/2*(a+b+(a-b)*exp(1j*pi*t))
    


###### Ipeng droit (attention ne comprend pas le dernier virage)

def ipengd(a,b,t):
    k=int(2*t)
    c=(b-1j*a)/(1-1j)
    L=[c,a,a+2*(b-a)*1j-(b-a)/2,b+(b-a)*sqrt(2)*exp(1j*pi/4),b]
    if t==1:return b #point d'arrivé
    else:
        tp=2*t-k
        l=(seg(a,c,tp),bezier(L,tp)) 
        return l[k]

#Ipeng gauche (attention ne comprend pas le dernier virage)
def ipengg(a,b,t):
    k=int(2*t)
    c=(a-1j*b)/(1-1j)
    L=[c,a,a+2*(a-b)*1j+(a-b)/2,b+(b-a)*sqrt(2)*exp(-1j*pi/4),b]
    if t==1:return b #point d'arrivé
    else:
        tp=2*t-k
        l=(seg(a,c,tp),bezier(L,tp)) 
        return l[k]
   

#itel droit : iteld
def iteld(a,b,t):
    z0=a+(1/2)*(1+1j)*(b-a)
    if t<1/3:return seg(a,b+1j*(b-a),3*t)
    if t>=1/3 and t<2/3:
        return exp(1j*(1-3*t)*3*pi/2)*(b+1j*(b-a)-z0)+z0
    else:
        return seg(a+1j*(b-a),b,3*t-2)

#itel gauche
def itelg(b,a,t):
    return iteld(a,b,1-t)

# Mli_basul droit de a à a+u

def mli_basuld(a,u,t):
    k=int(5*t)
    b=a+u
    if t==1:return a #point d'arrivé
    else:
        tp=5*t-k
        l=(vd(a,b,tp),pbd(b,-1j*u,tp),pbd(b,u,tp),pbd(b,1j*u,tp),vd(b,a,tp)) 
        return l[k]


#Mli basul droit de a à a+u
def mli_basulg(a,u,t):
        return mli_basuld(a,u,1-t)

### Rupture de direction à gauche et à droite (à gauche = on part à gauche)
def ruptureg(a,b,t):
    return a

def ruptured(a,b,t):
    return a

##### Carapace tortue #######
def carapaced(a,b,t): ##n= nombre de ligne à descendre
        u=(b-a)/abs(b-a)
        v=1j*u #repère local
        L1=[a,a-0.5*u+0.5*v,a-1*u,a-1*u-2.5*v,a+1*u-2.5*v]
        L2=list(reversed([b,b+0.5*u+0.5*v,b+1*u,b+1*u-2.5*v,b-1*u-2.5*v]))
        if t<=0.5:
                return bezier(L1,2*t)
        else:
                return bezier(L2,2*t-1)

def carapaceg(b,a,t):
        return carapaced(a,b,1-t)



def patted(a,b,t):
    L1=[a,a+0.3*sqrt(2)*(b-a)*exp(1j*pi/4),a+0.5*sqrt(2)*(b-a)*exp(1j*pi/4),\
        a+1.5*sqrt(2)*(b-a)*exp(1j*pi/4),b+(b-a)/2+(b-a)*sqrt(2)*exp(1j*pi/4)]
    L2=[b+(b-a)/2+(b-a)*sqrt(2)*exp(1j*pi/4),b+0.3*(b-a)*sqrt(2)*exp(1j*pi/4),b]
    if t<=0.5:
        return bezier(L1,2*t)
    else:
        return bezier(L2,2*t-1)

def pattedenvers(a,b,t):
    return patted(b,a,1-t)

def s(a,z):
        return a-(z-a).conjugate()

def patteg(a,b,t):
    return s(a,patted(a,b,t)) 



def pattegenvers(a,b,t):
    return patteg(b,a,1-t)

def iteld_souple(a,b,t):
    L1=[a,a+(1/4)*(b-a)*exp(1j*pi/4),b+(a-b)*(-1j),b+(1/4)*(b-a)*exp(1j*pi/4),b]
    L2=[a,a+(1/4)*(b-a)*exp(3j*pi/4),a+(b-a)*(1j),b+(1/4)*(b-a)*exp(3j*pi/4),b]
    if t<=1/3:
        return bezier(L1,3*t)
    elif t<=2/3:
        return vd(b, a, 3*t-1)
    else:
        return bezier(L2,3*t-2)

def itelg_souple(b,a,t):
        return iteld_souple(a,b,1-t)

def moins_id(a,u,t):
    k=1.5
    L=[a,a+0.5*u,a+k*u*exp(1j*pi/4),a+k*u*exp(-1j*pi/4),a+0.5*u,a]
    return bezier(L, t)

##### dictionnaire des valeurs de chaque mouvements

valeurs={seg:[1,1,1],
         vd:[1,0,-1j],
         vg:[0,1,1j],
         v_2_d:[2,0,-1j],
         v_2_g:[0,2,1j],
         Gvd:[0,-2,1j],
         Gvg:[-2,0,-1j],
         pbg:[0,0,-1j],
         pbd:[0,0,1j],
         pbg_concave:[0,0,-1j],
         pbd_concave:[0,0,1j],
         iteld:[1,0,-1j],
         itelg:[0,1,1j],
         iteld_souple:[1,0,-1j],
         itelg_souple:[0,1,1j],
         mli_basuld:[0,0,1j],
         mli_basulg:[0,0,1j],
         ipengd:[1,0,-1],
         ipengg:[0,1,-1],
         patted:[1,0,-1],
         patteg:[0,1,-1],
         pattedenvers:[-1,0,-1],
         pattegenvers:[0,-1,-1],
         vvd:[1,0,-1j],
         vvg:[0,1,1j],
         Gvvd:[2,0,-1j],
         Gvvg:[0,2,1j],
         carapaced:[0,-2,1j],
         carapaceg:[-2,0,-1j],
         ruptureg:[0,0,1j],
         ruptured:[0,0,-1j],
         moins_id:[0,0,-1]
         }

#### Calcul du couple [point,direction] suivant, c'est l'action de f sur (A,d)

def action(f,A,d):
    alpha,beta,gamma=valeurs[f][0],valeurs[f][1],valeurs[f][2]
    return [A+(alpha+1j*beta)*d,gamma*d]

#Création liste de deux listes :  points et directions
def actionliste(Lmouv,A0,d0):
    L,d=[A0],[d0]
    for m in Lmouv:
        nv_point=action(m,L[-1],d[-1])[0]
        nv_direction=action(m,L[-1],d[-1])[1]
        L+=[nv_point]
        d+=[nv_direction]
    return [L,d]


######Outils graphiques######


# Création de la courbe globale pour t in [0,1]

def f(lpoint,lmouv,ldirection,t): ## fonction auxiliaire
    L=[]
    l=len(lmouv)
    for k in range(l):
        if lmouv[k]==pbg or lmouv[k]==pbg_concave:
            f=lmouv[k]
            L+=[f(lpoint[k],ldirection[k]*exp(1j*pi/2),t)]
        elif lmouv[k]==pbd or lmouv[k]==pbd_concave:
            f=lmouv[k]
            L+=[f(lpoint[k],ldirection[k],t)]
        elif lmouv[k]==moins_id:
            f=lmouv[k]
            L+=[f(lpoint[k],ldirection[k]*exp(1j*pi/4),t)]
        else:
            L+=[lmouv[k](lpoint[k],lpoint[k+1],t)]
    return L

def chemin(lpoint,lmouv,ldirection,t):
    l=len(lmouv)
    k=int(l*t)
    if t==1:
        return lpoint[-1]
    else: 
        return f(lpoint,lmouv,ldirection,t*l-k)[k]



#Outil de grille : trace un segment (a,b)-(a1,b1)
def segment(a,b,a1,b1,color='k'):
    plt.plot([a,a1],[b,b1],color,linestyle='dashed')

def grille(l,c,pt_inf_gauche=0): #n lignes p colonnes
    a=pt_inf_gauche.real
    b=pt_inf_gauche.imag
    for i in range(0,c):
        segment(i+a,l-1+0.2+b,i+a,-0.2+b)
    for i in range(0,l):
        segment(-0.2+a,i+b,c-1+0.2+a,i+b)

#Test point du bord pour une grille n,p dont le coin inf gauche est donné
def bord(z,n,p,pt_inf_gauche):
    u,v=pt_inf_gauche.real,pt_inf_gauche.imag
    if z.real==u or (z.real==u+p-1) or (z.imag==v) or z.imag==v+n-1: 
        return True
    else:
        return False

#Test si un point est un coin
def coin(z,n,p,pt_inf_gauche):
    A0=pt_inf_gauche
    if z in {A0,A0+p-1,A0+(n-1)*1j,A0+p-1+(n-1)*1j}:
        return True
    else:
        return False




class tortue:
    #n,p,A0,d0 = taille de grille, point initial , direction départ
    #mouv_init = mouvement initial
    #A0 est un des coins de la grille
    #vertical dans le cas où la tortue est réalisée avec une rotation de 90°
    #Vertical reste à faire au 05/10/20
    def __init__(self,ligne=5,colonne=3,pt_inf_gauche=0,\
        A0=0,d0=1,mouv_init=seg,vertical=False):
        self.ligne=ligne
        self.colonne=colonne
        self.pt_inf_gauche=pt_inf_gauche
        self.A0=A0
        self.d0=d0
        self.mouv_init=mouv_init
        self.vertical=vertical
        #liste initial=e point de départ
        n=ligne
        p=colonne
        valeurs[Gvg][0]=-(p-1) # change selon la taille de grille
        valeurs[Gvd][1]=-(p-1) # change selon la taille de grille
        A1=action(mouv_init,A0,d0)[0]
        d1=action(mouv_init,A0,d0)[1]
        lpoint=[A0,A1]
        ldirection=[d0,d1]
        #liste initiale des déplacements
        ldepl=[mouv_init]
        while lpoint.count(A0)<3:
            Mk=lpoint[-1]
            dk=ldirection[-1]
            if not(bord(Mk,n,p,pt_inf_gauche)):
                ldepl.append(seg)
                lpoint.append(action(seg,Mk,dk)[0])
                ldirection.append(action(seg,Mk,dk)[1])

            elif coin(Mk,n,p,pt_inf_gauche):
                    if ldepl[-1]==seg:
                        if dk in {1,-1}:
                            ldepl.append(Gvg)
                            lpoint.append(action(Gvg,Mk,dk)[0])
                            ldirection.append(action(Gvg,Mk,dk)[1])
                        else:
                            ldepl.append(Gvd)
                            lpoint.append(action(Gvd,Mk,dk)[0])
                            ldirection.append(action(Gvd,Mk,dk)[1])
                    elif ldepl[-1] in {Gvg,Gvd}:
                        ldepl.append(seg)
                        lpoint.append(action(seg,Mk,dk)[0])
                        ldirection.append(action(seg,Mk,dk)[1])
                    elif ldepl[-1]==vd:
                        ldepl.append(vd)
                        lpoint.append(action(vd,Mk,dk)[0])
                        ldirection.append(action(vd,Mk,dk)[1])
                    else:
                        ldepl.append(vg)
                        lpoint.append(action(vg,Mk,dk)[0])
                        ldirection.append(action(vg,Mk,dk)[1])
            else: 
                if ldepl[-1]==seg:
                    if Mk.real==pt_inf_gauche.real or Mk.real==pt_inf_gauche.real+p-1:
                            if dk in {1,-1}: 
                                ldepl.append(vg)
                                lpoint.append(action(vg,Mk,dk)[0])
                                ldirection.append(action(vg,Mk,dk)[1])
                            else: 
                                ldepl.append(vd)
                                lpoint.append(action(vd,Mk,dk)[0])
                                ldirection.append(action(vd,Mk,dk)[1])
                    else:
                        if dk in {1,-1}: 
                                ldepl.append(vd)
                                lpoint.append(action(vd,Mk,dk)[0])
                                ldirection.append(action(vd,Mk,dk)[1])
                        else: 
                            ldepl.append(vg)
                            lpoint.append(action(vg,Mk,dk)[0])
                            ldirection.append(action(vg,Mk,dk)[1])
                            
                else:
                    ldepl.append(seg)
                    lpoint.append(action(seg,Mk,dk)[0])
                    ldirection.append(action(seg,Mk,dk)[1])
        self.lpoint=lpoint
        self.lmouv=ldepl 
        self.ldirection=ldirection
    
######## Tortue étoile ########


#Test si un point est un coin
def coin_haut(z,n,p,pt_inf_gauche):
    A0=pt_inf_gauche
    if z in {A0+(n-1)*1j,A0+p-1+(n-1)*1j}:
        return True
    else:
        return False

def coin_bas(z,n,p,pt_inf_gauche):
    A0=pt_inf_gauche
    if z in {A0,A0+p-1}:
        return True
    else:
        return False

def bas(z,n,p,pt_inf_gauche): ## En bas mais non coin
    A0=pt_inf_gauche
    if z.imag==A0.imag:
        return True
    else:
        return False

def haut(z,n,p,pt_inf_gauche): ## En haut
    A0=pt_inf_gauche
    if z.imag==A0.imag+n-1:
        return True
    else:
        return False

class tortue_etoile:
    #n,p,A0,d0 = taille de grille, point initial , direction départ
    #mouv_init = mouvement initial
    #A0 est un des coins de la grille
    #vertical dans le cas où la tortue est réalisée avec une rotation de 90°
    #Vertical reste à faire au 05/10/20
    def __init__(self,ligne=5,colonne=5,pt_inf_gauche=0+0j,\
        A0=0,d0=1,mouv_init=seg,vertical=False):
        self.ligne=ligne
        self.colonne=colonne
        self.pt_inf_gauche=pt_inf_gauche
        self.A0=A0
        self.d0=d0
        self.mouv_init=mouv_init
        self.vertical=vertical
        #liste initial=e point de départ
        n=ligne
        p=colonne
        valeurs[Gvg][0]=-(p-1) # change selon la taille de grille
        valeurs[Gvd][1]=-(p-1) # change selon la taille de grille
        A1=action(mouv_init,A0,d0)[0]
        d1=action(mouv_init,A0,d0)[1]
        lpoint=[A0,A1]
        ldirection=[d0,d1]
        #liste initiale des déplacements
        ldepl=[mouv_init]
        while lpoint.count(A0)<3:
            Mk=lpoint[-1]
            dk=ldirection[-1]
            if not(bord(Mk,n,p,pt_inf_gauche)):
                ldepl.append(seg)
            elif coin_haut(Mk,n,p,pt_inf_gauche):
                    if ldepl[-1]==seg:
                        if dk in {1,-1}:
                            ldepl.append(Gvg)
                        else:
                            ldepl.append(Gvd)
                    elif ldepl[-1] in {Gvg,Gvd}:
                        ldepl.append(seg)
                    elif ldepl[-1]==vd:
                        ldepl.append(vd)
                    else:
                        ldepl.append(vg)
            elif bas(Mk,n,p,pt_inf_gauche):
                if dk==-1:
                    ldepl.append(pbg)
                if dk==1:
                    if Mk.real==p-1:
                        ldepl.append(vg)
                    else:
                        ldepl.append(seg)
                if dk==1j:
                    if Mk.real==0:
                        ldepl.append(vd)
                    else:
                        ldepl.append(seg)
                if dk==-1j:
                    ldepl.append(pbd)
            else: 
                if ldepl[-1]==seg:
                    if Mk.real==pt_inf_gauche.real or Mk.real==pt_inf_gauche.real+p-1:
                            if dk in {1,-1}: 
                                ldepl.append(vg)
                            else: 
                                ldepl.append(vd)
                    else:
                        if dk in {1,-1}: 
                                ldepl.append(vd)
                        else: 
                            ldepl.append(vg)     
                else:
                    ldepl.append(seg)
            lpoint.append(action(ldepl[-1],Mk,dk)[0])
            ldirection.append(action(ldepl[-1],Mk,dk)[1])
        self.lpoint=[M+pt_inf_gauche for M in lpoint]
        self.lmouv=ldepl 
        self.ldirection=ldirection


class billard:
    #n,p,A0,d0 = taille de grille, point initial , direction départ
    #mouv_init = mouvement initial
    #A0 est un des coins de la grille
    #vertical dans le cas où la tortue est réalisée avec une rotation de 90°
    #Vertical reste à faire au 05/10/20
    def __init__(self,ligne=5,colonne=5,pt_inf_gauche=0+0j,\
        A0=0,d0=1,mouv_init=seg,vertical=False):
        self.ligne=ligne
        self.colonne=colonne
        self.pt_inf_gauche=pt_inf_gauche
        self.A0=A0
        self.d0=d0
        self.mouv_init=mouv_init
        self.vertical=vertical
        #liste initial=e point de départ
        n=ligne
        p=colonne
        valeurs[Gvg][0]=-(p-1) # change selon la taille de grille
        valeurs[Gvd][1]=-(p-1) # change selon la taille de grille
        A1=action(mouv_init,A0,d0)[0]
        d1=action(mouv_init,A0,d0)[1]
        lpoint=[A0,A1]
        ldirection=[d0,d1]
        #liste initiale des déplacements
        ldepl=[mouv_init]
        while lpoint.count(A0)<3:
            Mk=lpoint[-1]
            dk=ldirection[-1]
            if not(bord(Mk,n,p,pt_inf_gauche)):
                ldepl.append(seg)
            elif bas(Mk,n,p,pt_inf_gauche):
                if dk==-1:
                    ldepl.append(pbg)
                if dk==1:
                    if Mk.real==p-1:
                        ldepl.append(vg)
                    else:
                        ldepl.append(seg)
                if dk==1j:
                    if Mk.real==0:
                        ldepl.append(vd)
                    else:
                        ldepl.append(seg)
                if dk==-1j:
                    ldepl.append(pbd)
            elif haut(Mk,n,p,pt_inf_gauche):
                if dk==-1:
                    if Mk.real==0:
                        ldepl.append(vg)
                    else:
                        ldepl.append(seg)
                elif dk==1:
                    ldepl.append(pbg)
                elif dk==1j:
                    ldepl.append(pbd)
                else:
                    if Mk.real==p-1:
                        ldepl.append(vd)
                    else:
                        ldepl.append(seg)
            else: 
                if ldepl[-1]==seg:
                    if Mk.real==0 or Mk.real==p-1:
                            if dk in {1,-1}: 
                                ldepl.append(vg)
                            else: 
                                ldepl.append(vd)
                    else:
                        if dk in {1,-1}: 
                                ldepl.append(vd)
                        else: 
                            ldepl.append(vg)     
                else:
                    ldepl.append(seg)
            lpoint.append(action(ldepl[-1],Mk,dk)[0])
            ldirection.append(action(ldepl[-1],Mk,dk)[1])
        self.lpoint=[M+pt_inf_gauche for M in lpoint]
        self.lmouv=ldepl 
        self.ldirection=ldirection


class tortue_etoile_twist:
    #n,p,A0,d0 = taille de grille, point initial , direction départ
    #mouv_init = mouvement initial
    #A0 est un des coins de la grille
    # seuil et seuil+1 : là où l'on fait les twists
    def __init__(self,seuil=3,ligne=5,colonne=5,pt_inf_gauche=0+0j,\
        A0=0,d0=1,mouv_init=seg):
        self.seuil=seuil
        self.ligne=ligne
        self.colonne=colonne
        self.pt_inf_gauche=pt_inf_gauche
        self.A0=A0
        self.d0=d0
        self.mouv_init=mouv_init
        #liste initial=e point de départ
        n=ligne
        p=colonne
        valeurs[Gvg][0]=-(p-1) # change selon la taille de grille
        valeurs[Gvd][1]=-(p-1) # change selon la taille de grille
        A1=action(mouv_init,A0,d0)[0]
        d1=action(mouv_init,A0,d0)[1]
        lpoint=[A0,A1]
        ldirection=[d0,d1]
        #liste initiale des déplacements
        ldepl=[mouv_init]
        while lpoint.count(A0)<3:
            Mk=lpoint[-1]
            dk=ldirection[-1]
            if Mk.imag==seuil and dk==1j:
                    ldepl+=[vd]
            elif Mk.imag==seuil and dk==1:
                    ldepl+=[vg]
            elif Mk.imag==seuil+1 and dk==-1:
                    ldepl+=[vg]
            elif Mk.imag==seuil+1 and dk==-1j:
                    ldepl+=[vd]
            elif not(bord(Mk,n,p,pt_inf_gauche)):
                ldepl.append(seg)
            elif coin_haut(Mk,n,p,pt_inf_gauche):
                    if ldepl[-1]==seg:
                        if dk in {1,-1}:
                            ldepl.append(Gvg)
                        else:
                            ldepl.append(Gvd)
                    elif ldepl[-1] in {Gvg,Gvd}:
                        ldepl.append(seg)
                    elif ldepl[-1]==vd:
                        ldepl.append(vd)
                    else:
                        ldepl.append(vg)
            elif bas(Mk,n,p,pt_inf_gauche):
                if dk==-1:
                    ldepl.append(pbg)
                if dk==1:
                    if Mk.real==p-1:
                        ldepl.append(vg)
                    else:
                        ldepl.append(seg)
                if dk==1j:
                    if Mk.real==0:
                        ldepl.append(vd)
                    else:
                        ldepl.append(seg)
                if dk==-1j:
                    ldepl.append(pbd)
            elif ldepl[-1]==vg and Mk.real==pt_inf_gauche.real and dk==1j:
                ldepl+=[vd]
            elif ldepl[-1]==vg and Mk.real==pt_inf_gauche.real+p-1 and dk==-1j:
                ldepl+=[vd]
            elif ldepl[-1]==vd and Mk.real==pt_inf_gauche.real and dk==-1:
                ldepl+=[vg]
            elif ldepl[-1]==vd and Mk.real==pt_inf_gauche.real+p-1 and dk==1:
                ldepl+=[vg]
            else: 
                if ldepl[-1]==seg:
                    if Mk.real==pt_inf_gauche.real or Mk.real==pt_inf_gauche.real+p-1:
                            if dk in {1,-1}: 
                                ldepl.append(vg)
                            else: 
                                ldepl.append(vd)
                    else:
                        if dk in {1,-1}: 
                                ldepl.append(vd)
                        else: 
                            ldepl.append(vg)     
                elif Mk.imag==seuil and Mk.real==0 and dk==-1:
                    ldepl+=[vg]
                elif Mk.imag==seuil and Mk.real==0 and dk==-1j:
                    ldepl+=[seg]
                elif Mk.imag==seuil and Mk.real==p-1 and dk==-1j:
                    ldepl+=[vd]
                elif Mk.imag==seuil+1 and Mk.real==p-1 and dk==1j:
                    ldepl+=[seg]
                elif Mk.imag==seuil+1 and Mk.real==p-1 and dk==1:
                    ldepl+=[vg]
                elif Mk.imag==seuil+1 and Mk.real==0 and dk==1j:
                    ldepl+=[vd]
                else:
                    ldepl.append(seg)
            lpoint.append(action(ldepl[-1],Mk,dk)[0])
            ldirection.append(action(ldepl[-1],Mk,dk)[1])
        self.lpoint=[M+pt_inf_gauche for M in lpoint]
        self.lmouv=ldepl 
        self.ldirection=ldirection

class billard_twist:
    #n,p,A0,d0 = taille de grille, point initial , direction départ
    #mouv_init = mouvement initial
    #A0 est un des coins de la grille
    #Reprend l'algo B mais place des miroirs verticaux entre la 
    #ligne seuil et seuil+1
    def __init__(self,ligne=5,colonne=5,seuil=3,pt_inf_gauche=0+0j,\
        A0=0,d0=1,mouv_init=seg):
        self.ligne=ligne
        self.colonne=colonne
        self.pt_inf_gauche=pt_inf_gauche
        self.A0=A0
        self.d0=d0
        self.mouv_init=mouv_init
        #liste initial=e point de départ
        n=ligne
        p=colonne
        valeurs[Gvg][0]=-(p-1) # change selon la taille de grille
        valeurs[Gvd][1]=-(p-1) # change selon la taille de grille
        A1=action(mouv_init,A0,d0)[0]
        d1=action(mouv_init,A0,d0)[1]
        lpoint=[A0,A1]
        ldirection=[d0,d1]
        #liste initiale des déplacements
        ldepl=[mouv_init]
        while lpoint.count(A0)<3:
            Mk=lpoint[-1]
            dk=ldirection[-1]
            if Mk.imag==seuil and dk==1j:
                    ldepl+=[vd]
            elif Mk.imag==seuil and dk==1:
                    ldepl+=[vg]      
            elif Mk.imag==seuil+1 and dk==-1:
                    ldepl+=[vg]
            elif Mk.imag==seuil+1 and dk==-1j:
                    ldepl+=[vd]
            elif not(bord(Mk,n,p,pt_inf_gauche)):
                ldepl.append(seg)
            elif bas(Mk,n,p,pt_inf_gauche):
                if dk==-1:
                    ldepl.append(pbg)
                if dk==1:
                    if Mk.real==p-1:
                        ldepl.append(vg)
                    else:
                        ldepl.append(seg)
                if dk==1j:
                    if Mk.real==0:
                        ldepl.append(vd)
                    else:
                        ldepl.append(seg)
                if dk==-1j:
                    ldepl.append(pbd)
            elif haut(Mk,n,p,pt_inf_gauche):
                if dk==-1:
                    if Mk.real==0:
                        ldepl.append(vg)
                    else:
                        ldepl.append(seg)
                elif dk==1:
                    ldepl.append(pbg)
                elif dk==1j:
                    ldepl.append(pbd)
                else:
                    if Mk.real==p-1:
                        ldepl.append(vd)
                    else:
                        ldepl.append(seg)
            else: 
                if ldepl[-1]==seg:
                    if Mk.real==0 or Mk.real==p-1:
                            if dk in {1,-1}: 
                                ldepl.append(vg)
                            else: 
                                ldepl.append(vd)
                    else:
                        if dk in {1,-1}: 
                                ldepl.append(vd)
                        else: 
                            ldepl.append(vg)
                elif Mk.imag==seuil and Mk.real==0 and dk==-1:
                    ldepl+=[vg]
                elif Mk.imag==seuil and Mk.real==0 and dk==-1j:
                    ldepl+=[seg]
                elif Mk.imag==seuil and Mk.real==p-1 and dk==-1j:
                    ldepl+=[vd]
                elif Mk.imag==seuil+1 and Mk.real==p-1 and dk==1j:
                    ldepl+=[seg]
                elif Mk.imag==seuil+1 and Mk.real==p-1 and dk==1:
                    ldepl+=[vg]
                elif Mk.imag==seuil+1 and Mk.real==0 and dk==1j:
                    ldepl+=[vd]
                else:
                    ldepl.append(seg)
            lpoint.append(action(ldepl[-1],Mk,dk)[0])
            ldirection.append(action(ldepl[-1],Mk,dk)[1])
        self.lpoint=[M+pt_inf_gauche for M in lpoint]
        self.lmouv=ldepl 
        self.ldirection=ldirection

class billard_twist_horizontal:
    #n,p,A0,d0 = taille de grille, point initial , direction départ
    #mouv_init = mouvement initial
    #A0 est un des coins de la grille
    #Reprend l'algo B mais place des miroirs horizontaux entre la 
    #colonne seuil et seuil+1
    def __init__(self,ligne=5,colonne=5,seuil=3,pt_inf_gauche=0+0j,\
        A0=0,d0=1,mouv_init=seg):
        self.ligne=ligne
        self.colonne=colonne
        self.pt_inf_gauche=pt_inf_gauche
        self.A0=A0
        self.d0=d0
        self.mouv_init=mouv_init
        #liste initial=e point de départ
        n,l=ligne,ligne
        p,c=colonne,colonne
        valeurs[Gvg][0]=-(p-1) # change selon la taille de grille
        valeurs[Gvd][1]=-(p-1) # change selon la taille de grille
        A1=action(mouv_init,A0,d0)[0]
        d1=action(mouv_init,A0,d0)[1]
        lpoint=[A0,A1]
        ldirection=[d0,d1]
        #liste initiale des déplacements
        ldepl=[mouv_init]
        while lpoint.count(A0)<3:
            Mk=lpoint[-1]
            dk=ldirection[-1]
            if (Mk.real==seuil and dk==1):
                if Mk.imag==l-1:
                    ldepl+=[pbg]
                else:
                    ldepl+=[vd]
            elif Mk.real==seuil and dk==-1j:
                if Mk.imag==0:
                    ldepl+=[pbd]
                else:
                    ldepl+=[vg]      
            elif Mk.real==seuil+1 and dk==1j:
                if Mk.imag==l-1:
                    ldepl+=[pbd]
                else:
                    ldepl+=[vg]
            elif Mk.real==seuil+1 and dk==-1:
                if Mk.imag==0:
                    ldepl+=[pbg]
                else:
                    ldepl+=[vd]
            elif not(bord(Mk,n,p,pt_inf_gauche)):
                ldepl.append(seg)
            elif bas(Mk,n,p,pt_inf_gauche):
                if dk==-1:
                    ldepl.append(pbg)
                if dk==1:
                    if Mk.real==p-1:
                        ldepl.append(vg)
                    else:
                        ldepl.append(seg)
                if dk==1j:
                    if Mk.real==0:
                        ldepl.append(vd)
                    else:
                        ldepl.append(seg)
                if dk==-1j:
                    ldepl.append(pbd)
            elif haut(Mk,n,p,pt_inf_gauche):
                if dk==-1:
                    if Mk.real==0:
                        ldepl.append(vg)
                    else:
                        ldepl.append(seg)
                elif dk==1:
                    ldepl.append(pbg)
                elif dk==1j:
                    ldepl.append(pbd)
                else:
                    if Mk.real==p-1:
                        ldepl.append(vd)
                    else:
                        ldepl.append(seg)
            else: 
                if ldepl[-1]==seg:
                    if Mk.real==0 or Mk.real==p-1:
                            if dk in {1,-1}: 
                                ldepl.append(vg)
                            else: 
                                ldepl.append(vd)
                    else:
                        if dk in {1,-1}: 
                                ldepl.append(vd)
                        else: 
                            ldepl.append(vg)
                else:
                    ldepl.append(seg)
            lpoint.append(action(ldepl[-1],Mk,dk)[0])
            ldirection.append(action(ldepl[-1],Mk,dk)[1])
        self.lpoint=[M+pt_inf_gauche for M in lpoint]
        self.lmouv=ldepl 
        self.ldirection=ldirection


def dessin(lmouv,A0,d0,ligne=3,colonne=3,N=3000,epaisseur=3,couleur='b',\
    A0_visible=False,d0_visible=False):
    l,c=ligne,colonne
    lpoint=actionliste(lmouv,A0,d0)[0]
    ldirection=actionliste(lmouv,A0,d0)[1]
    sub=numpy.linspace(0,1,N)
    x=[chemin(lpoint,lmouv,ldirection,t).real for t in sub]
    y=[chemin(lpoint,lmouv,ldirection,t).imag for t in sub]
    fig, ax = plt.subplots() #fig et ax sont des classes
    line,=ax.plot(x,y,linewidth=3,color='b',path_effects=[path_effects.SimpleLineShadow(),\
                       path_effects.Normal(offset=1)])
    grille(l,c)
    if A0_visible==True:
        plt.scatter([A0.real], [A0.imag],c = ['lawgngreen'],s = [110],marker = 'o')
    ax.grid(False) #quadrillage
    ax.axis([-2,c+2,-2,l+2]) #on fixe les xmin,xmax,ymin,ymax
    ax.set_axis_off()
    ax.set_xticklabels([]) #pas de label en x
    ax.set_yticklabels([]) #pas de label en y
    ax.tick_params(axis='x', colors=(0,0,0,0)) #pas de marqueur x
    ax.tick_params(axis='y', colors=(0,0,0,0)) #pas de marqueur y
    plt.show()

