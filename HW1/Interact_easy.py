import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from sympy import symbols, Eq, solve
from Structure import Triangle, Voronoi_cell
from scipy.spatial import Delaunay # à enlever après

#######################################################################################################################
#DELANNAY DEJA TOUT FAIT 
#########################################################################################################################""

# --- Initialisation de la figure ---
fig, ax = plt.subplots()
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_title("Cliquez pour ajouter des points (animation + Delaunay)")

# --- Données ---
points_x, points_y = [1,1,2,2,3,2], [1,2,1,2,2,3]  # Listes des coordonnées des points

# --- Tracés initiaux ---
(scatter,) = ax.plot([], [], "ro")
tri_lines = None  # objet pour la triangulation

# --- Fonction de mise à jour de l'animation ---
def update(frame):
    global tri_lines
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    arr_x = np.array(points_x, dtype=np.double)
    arr_y = np.array(points_y, dtype=np.double)

    # Nettoyer la triangulation précédente
    if tri_lines:
        for t in tri_lines:
            t.remove()

    # Mettre à jour les points
    scatter.set_data(points_x, points_y)

    # Dessiner la nouvelle triangulation si assez de points
    if len(points_x) >= 3:
        tri = Delaunay(np.array(list(zip(arr_x, arr_y))))
        tri_lines = ax.triplot(arr_x, arr_y, tri.simplices, "b-", lw=1)
    else:
        tri_lines = []

    return scatter,

# --- Fonction de clic ---
def onclick(event):
    if event.inaxes != ax:
        return
    # Ajout d’un point
    points_x.append(event.xdata)
    points_y.append(event.ydata)
    print(f"Point ajouté : ({event.xdata:.2f}, {event.ydata:.2f})")

# --- Connexion du clic ---
fig.canvas.mpl_connect("button_press_event", onclick)

# --- Animation ---
ani = animation.FuncAnimation(fig, update, interval=100, blit=False)

plt.show()


############################################################################################
#Voronoi Cells : Dual Graph of the Delaunay Triangulation
#The vertices of the Voronoi diagram are the circumcenters of the triangles in the Delaunay triangulation
#The midpoints of the Voronoi cells are the vertices of the Delaunay triangulation
#############################################################################################


def circumcenter(a, b, c):
        x, y , R = symbols('x y R')
        eq1 = Eq((a[0]-x)**2 + (a[1]-y)**2, R**2)
        eq2 = Eq((b[0]-x)**2 + (b[1]-y)**2, R**2)
        eq3 = Eq((c[0]-x)**2 + (c[1]-y)**2, R**2)

        sol = solve((eq1, eq2, eq3), (x, y, R))
        return (float(sol[0][0]), float(sol[0][1]))


def Voronoi(pts): #crée une liste de la structure Voronoi_cells dans Structure
    
    i = 0
    triangles = []
    for trig in Delaunay(np.array(pts)).simplices:    #à changer ICI
        tri = Triangle(pts[trig[0]], pts[trig[1]], pts[trig[2]], i)
        i += 1
        tri.center = circumcenter(pts[trig[0]], pts[trig[1]], pts[trig[2]])   #Chaque coin d'une cell est le centre du cercle circonscrit du triangle de Delannay
        triangles.append(tri)
    
    Voronoi_cells = []

    for p in pts:
        cell = Voronoi_cell(p, [])
        corners = []
        for tri in triangles:
            if p in tri.points : 
                corners = append_unique(corners, tri.center) #case where two triangles share the same circumcenter
                #print(f"Point {p} is in Triangle {tri.j} with center at {tri.center}")
        corners = sorted(corners, key=lambda point: np.arctan2(point[1]-p[1], point[0]-p[0]))
        cell.corners = corners
        Voronoi_cells.append(cell)
    return Voronoi_cells

#checks if a center is already in corners (with a tolerance) [case where two triangles share the same circumcenter]
def append_unique(corners, center, tol=1e-9):
    for c in corners:
        if np.linalg.norm(np.array(c) - np.array(center)) < tol: #norme euclidienne
            return corners
    corners.append(center)
    return corners

def prod_vectoriel (a, b, c):
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])


#on renvoit une liste de tout les centres des cells qui appartiennent au convex Hull
#Gift Wrapping (Algorithme de Jarvis March)
def convex_hull(pts):  
    hull=[]

    l = pts[0]
    p=l
    while True:
        hull.append(p)
        q = pts[0]
        for r in pts:
            if (q == p) or (prod_vectoriel (p, q, r) < 0):
                q = r
        p = q
        if p == l:
            break
    return hull



# --- Initialisation de la figure Voronoi---
fig, ax = plt.subplots()
ax.set_aspect("equal") 
xmin, xmax = 0, 10
ymin, ymax = 0, 10
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_title("Cliquez pour ajouter des points (Voronoi Cells)")

# --- Données ---
pts = [ (5,3), (2,7), (8,4), (6,6), (3,3), (6,7), (4, 5) ]  # Initial test points
pts = [ (2,9), (3,8), (4,5)]  # Initial test points
pts = sorted(pts, key=lambda point: (point[0], point[1]))



# Plot setup
ax.set_title("Voronoi Diagram (points test)")

def isincell(c, A, B): #check si un point fait partie d'une des 2 cells ayant pour centre A ou B
        for cl in cells:
            if c in cl.corners and A == cl.center:
                return True
            if c in cl.corners and B == cl.center:
                return True
        return False

def point_on_line(c, M, n, tol=1e-6): #return True si le point c fait partie de la droite selon le point M et le vecteur n
        return abs((c[0] - M[0]) * n[1] - (c[1] - M[1]) * n[0]) < tol

# Plot the Voronoi cells
if len(pts) >= 3: #>=3 forme au moins un triangle de Delannay

    cells = Voronoi(pts)
    hull_points = convex_hull(pts)
    print(hull_points)
    for cell in cells:
        corners = np.array(cell.corners)
        if cell.center in hull_points and len(corners) > 2: #ne dessine que les points appartenant au convex Hull (= segments infinies)
            ax.plot(cell.center[0], cell.center[1], 'ro')  # cell center
            ax.plot(corners[:, 0], corners[:, 1], 'go')  # cell corners
            continue
        #dessine le reste
        ax.plot(corners[:, 0], corners[:, 1], 'b-')  # cell boundary
        ax.plot(cell.center[0], cell.center[1], 'ro')  # cell center
        ax.plot(corners[:, 0], corners[:, 1], 'go')  # cell corners
            
    #infinite segments
    center = np.mean(pts, axis=0)
    data=[]

    #Par sureté : regarder si le convex Hull est dans le sens anti-horaire ou horaire = > ça change l'orientation du vecteur normal n après
    signed_area = 0
    for i in range(len(hull_points)):
        A = hull_points[i]
        B = hull_points[(i+1) % len(hull_points)]
        signed_area += (B[0] - A[0]) * (B[1] + A[1])  #somme des produits vectoriels ? un truc du genre
    ccw = signed_area < 0  # True = horaire, False = anti-horaire

    for i in range(len(hull_points)):
        A = hull_points[i]
        B = hull_points[(i+1) % len(hull_points)]
        C = hull_points[(i+2) % len(hull_points)]

        # Determiner le milieu M et le vecteur directeur v
        M = np.array([(A[0] + B[0])/2, (A[1] + B[1])/2])
        v = np.array([B[0] - A[0], B[1] - A[1]])
        v = v / np.linalg.norm(v)

        if abs(prod_vectoriel(A, B, C)) < 1e-9: #vérifie si les points sont alignés (on sait jamais)
            continue

        # Vecteur Normale au segment AB
        n = np.array([-v[1], v[0]]) if ccw else np.array([v[1], -v[0]])

        # Orientation vers l’extérieur
        if np.dot(M - center, n) < 0:
            n = -n

        # Tracer la demi-droite infinie
        
        data.append((M, n, A, B))
    

    print("Data for infinite edges :", data)
    
    cs = []
    for c in cells:
        for corner in c.corners:
            cs.append(corner)
     #print(cs) #corners

    for M, n, A, B in data:
        for c in cs:
            if len(cs)==1: #imagine il a qu'un bord 
                far = (c[0] + n[0] * 12, c[1] + n[1] * 12)  # point lointain pour simuler l'infini
                plt.plot([c[0], far[0]], [c[1], far[1]], 'g--')
                break
            if point_on_line(c, M, n)  and isincell(c, A, B):
                print("Intersection found at corner :", c)
                far = (c[0] + n[0] * 12, c[1] + n[1] * 12)  # point lointain pour simuler l'infini
                plt.plot([c[0], far[0]], [c[1], far[1]], 'g--')
            #print("Infinite edge from", M, "in direction", n)
            


if len(pts)<3:
    a = pts[0]
    b = pts[1]
    M = np.array([(a[0] + b[0])/2, (a[1] + b[1])/2])
    v = np.array([b[0] - a[0], b[1] - a[1]])
    v = v / np.linalg.norm(v)
    n = np.array([-v[1], v[0]])
    ax.scatter([p[0] for p in pts], [p[1] for p in pts], color='red')  # cell centers
    plt.plot([M[0]-n[0] * 12, M[0] + n[0] * 12], [M[1]-n[1] * 12, M[1] + n[1] * 12], 'g--')




# Plot original points
#ax.scatter([p[0] for p in pts], [p[1] for p in pts], color='black', label='Input points')

ax.legend()
plt.show()

