#!/usr/bin/env python3
import MDAnalysis as mda
import numpy as np
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # niezbędne do wykresów 3D
import os
from matplotlib import cm
import random

# Ustawienia parametrów analizy
eps_distance = 3.0       # Zmniejszony promień sąsiedztwa w DBSCAN [Å]
min_cluster_size = 100   # Minimalny rozmiar interesującego klastra
max_cluster_size = 5000  # Maksymalny rozmiar klastra (aby uniknąć całej wody jako jeden klaster)
min_samples = 10         # Zwiększona liczba punktów do zdefiniowania klastra
max_clusters_to_show = 3  # Liczba klastrów do pokazania
max_atoms_per_cluster = 500  # Maksymalna liczba atomów do wyświetlenia z każdego klastra

# Ścieżki do plików w katalogu wyżej (tip4p_273K)
# Używamy ścieżki względnej od miejsca uruchomienia skryptu
tpr_file = os.path.join("..", "..", "tip4p_273K", "md.tpr")
xtc_file = os.path.join("..", "..", "tip4p_273K", "md.xtc")

# Sprawdzenie, czy pliki istnieją
if not os.path.exists(tpr_file):
    print(f"Błąd: Plik {tpr_file} nie istnieje!")
    print(f"Bieżący katalog: {os.getcwd()}")
    exit(1)
if not os.path.exists(xtc_file):
    print(f"Błąd: Plik {xtc_file} nie istnieje!")
    exit(1)

# Wczytanie trajektorii
u = mda.Universe(tpr_file, xtc_file)

# Selektor atomów wody – sprawdź czy nazwa reszt jest zgodna z Twoimi danymi
water = u.select_atoms("resname SOL or resname TIP4")

print("Liczba atomów wody:", len(water))

# Utworzenie katalogu na wyniki, jeśli nie istnieje
output_dir = "cluster_results"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Funkcja do wizualizacji tylko interesujących klastrów
def visualize_interesting_clusters(coords, labels, frame_num, interesting_clusters, cluster_sizes):
    """
    Wizualizacja tylko interesujących klastrów z ograniczoną liczbą atomów.
    
    Parameters:
    -----------
    coords : numpy.ndarray
        Współrzędne atomów
    labels : numpy.ndarray
        Etykiety klastrów dla każdego atomu
    frame_num : int
        Numer klatki
    interesting_clusters : list
        Lista etykiet interesujących klastrów
    cluster_sizes : dict
        Słownik z rozmiarami klastrów
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Pokaż tylko wybrane interesujące klastry
    for i, label in enumerate(interesting_clusters):
        # Wybierz punkty należące do danego klastra
        mask = labels == label
        cluster_points = coords[mask]
        cluster_size = cluster_sizes[label]
        
        # Jeśli klaster jest duży, wybierz losową próbkę punktów
        if len(cluster_points) > max_atoms_per_cluster:
            indices = random.sample(range(len(cluster_points)), max_atoms_per_cluster)
            cluster_points = cluster_points[indices]
        
        # Użyj różnych kolorów dla każdego klastra
        color = cm.tab10(i % 10)
        
        # Dodaj klaster do wykresu
        ax.scatter(
            cluster_points[:, 0], 
            cluster_points[:, 1], 
            cluster_points[:, 2],
            color=color,
            s=15,  # Większy rozmiar punktów dla lepszej widoczności
            alpha=0.8,
            marker='o',
            label=f'Klaster {label} ({cluster_size} atomów)'
        )
    
    # Dodaj legendę
    ax.legend()
    
    # Dodaj tytuł i etykiety osi
    plt.title(f"Frame {frame_num} - Interesujące formacje wody")
    ax.set_xlabel("X [Å]")
    ax.set_ylabel("Y [Å]")
    ax.set_zlabel("Z [Å]")
    
    # Równe skalowanie osi dla lepszej perspektywy
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
    
    # Ustaw równe skalowanie
    max_range = 0.5 * max([x_range, y_range, z_range])
    ax.set_xlim3d([x_middle - max_range, x_middle + max_range])
    ax.set_ylim3d([y_middle - max_range, y_middle + max_range])
    ax.set_zlim3d([z_middle - max_range, z_middle + max_range])
    
    # Zapisz wizualizację do pliku
    plt.savefig(os.path.join(output_dir, f"interesting_clusters_frame_{frame_num}.png"), dpi=300)
    
    return fig

# Funkcja do testowania różnych parametrów DBSCAN
def find_optimal_parameters(coords, start_eps=2.0, end_eps=4.0, step=0.5):
    """
    Testuje różne wartości eps dla DBSCAN i zwraca najlepsze parametry.
    """
    best_eps = None
    best_n_clusters = 0
    
    for eps in np.arange(start_eps, end_eps + step, step):
        clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
        labels = clustering.labels_
        unique_labels = set(labels) - {-1}
        n_clusters = len(unique_labels)
        
        # Oblicz rozmiary klastrów
        cluster_sizes = {}
        for label in unique_labels:
            cluster_sizes[label] = np.sum(labels == label)
        
        # Sprawdź, czy znaleziono interesujące klastry (nie za duże, nie za małe)
        interesting_clusters = 0
        for size in cluster_sizes.values():
            if min_cluster_size <= size <= max_cluster_size:
                interesting_clusters += 1
        
        print(f"  Testowanie eps={eps:.1f}: znaleziono {n_clusters} klastrów, {interesting_clusters} interesujących")
        
        # Jeśli znaleziono więcej interesujących klastrów niż poprzednio, zaktualizuj najlepsze parametry
        if interesting_clusters > best_n_clusters:
            best_n_clusters = interesting_clusters
            best_eps = eps
    
    return best_eps if best_eps is not None else eps_distance

# Iteracja po wszystkich klatkach trajektorii
for ts in u.trajectory:
    # Pobranie współrzędnych atomów w bieżącej klatce (macierz N x 3)
    coords = water.positions
    
    # Raportuj co 100 klatek
    if ts.frame % 100 == 0:
        print(f"Analizowanie klatki {ts.frame}...")
    
    # Najpierw spróbuj z domyślnymi parametrami
    clustering = DBSCAN(eps=eps_distance, min_samples=min_samples).fit(coords)
    labels = clustering.labels_
    
    # Wyznaczenie unikalnych klastrów (pomijamy szum - etykieta -1)
    unique_labels = set(labels) - {-1}
    n_clusters = len(unique_labels)
    
    # Obliczanie rozmiarów poszczególnych klastrów
    cluster_sizes = {}
    for label in unique_labels:
        cluster_sizes[label] = np.sum(labels == label)
    
    # Sortowanie klastrów według rozmiaru (od największego)
    sorted_clusters = sorted(cluster_sizes.items(), key=lambda x: x[1], reverse=True)
    
    # Pobierz największy rozmiar klastra
    max_cluster_size_found = sorted_clusters[0][1] if sorted_clusters else 0
    
    # Jeśli znaleziono tylko jeden duży klaster (cała woda), spróbuj dostosować parametry
    if n_clusters == 1 and max_cluster_size_found > max_cluster_size:
        if ts.frame % 100 == 0:
            print(f"  Znaleziono jeden duży klaster ({max_cluster_size_found} atomów). Próba dostosowania parametrów...")
        
        # Spróbuj znaleźć lepsze parametry dla tej klatki
        optimal_eps = find_optimal_parameters(coords)
        
        # Uruchom DBSCAN ponownie z lepszymi parametrami
        clustering = DBSCAN(eps=optimal_eps, min_samples=min_samples).fit(coords)
        labels = clustering.labels_
        
        # Zaktualizuj informacje o klastrach
        unique_labels = set(labels) - {-1}
        n_clusters = len(unique_labels)
        
        cluster_sizes = {}
        for label in unique_labels:
            cluster_sizes[label] = np.sum(labels == label)
        
        sorted_clusters = sorted(cluster_sizes.items(), key=lambda x: x[1], reverse=True)
        max_cluster_size_found = sorted_clusters[0][1] if sorted_clusters else 0
        
        if ts.frame % 100 == 0:
            print(f"  Po dostosowaniu: eps={optimal_eps:.1f}, liczba klastrów = {n_clusters}, największy = {max_cluster_size_found}")
    
    # Raportuj co 100 klatek
    if ts.frame % 100 == 0:
        print(f"Frame = {ts.frame}, liczba klastrów = {n_clusters}, największy klaster = {max_cluster_size_found} atomów")
    
    # Znajdź interesujące klastry (o odpowiednim rozmiarze)
    interesting_clusters = []
    for label, size in sorted_clusters:
        if min_cluster_size <= size <= max_cluster_size:
            interesting_clusters.append(label)
            if len(interesting_clusters) >= max_clusters_to_show:
                break
    
    # Jeśli znaleziono interesujące klastry, zapisz klatkę i wizualizuj
    if interesting_clusters:
        print(f"UWAGA! Znaleziono interesujące klastry w frame = {ts.frame}")
        for label in interesting_clusters:
            print(f"  - Klaster {label}: {cluster_sizes[label]} atomów")
        
        # Zapisz wszystkie atomy wody w tej klatce do pliku PDB
        output_file = os.path.join(output_dir, f"interesting_clusters_frame_{ts.frame}.pdb")
        water.write(output_file)
        
        # Wizualizuj interesujące klastry
        fig = visualize_interesting_clusters(coords, labels, ts.frame, interesting_clusters, cluster_sizes)
        plt.show()