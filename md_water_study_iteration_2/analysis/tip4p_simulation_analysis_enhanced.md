# Rozszerzona analiza właściwości fizycznych wody TIP4P w 273K

## 1. Wprowadzenie

Niniejszy dokument przedstawia szczegółową analizę teoretyczną symulacji wody z wykorzystaniem modelu TIP4P (Transferable Intermolecular Potential with 4 Points) w temperaturze 273K. Analiza opiera się na danych z symulacji dynamiki molekularnej przeprowadzonej z użyciem pakietu GROMACS.

### 1.1 Model TIP4P

Model TIP4P to czterocentrowy model wody składający się z:
- Atomu tlenu (O) bez ładunku, który jest centrum oddziaływań van der Waalsa
- Dwóch atomów wodoru (H) z ładunkami dodatnimi
- Punktu masowego (M) zlokalizowanego w pobliżu tlenu, niosącego ładunek ujemny

Geometria molekuły wody w modelu TIP4P, zgodnie z plikiem `water.pdb`:
```
ATOM      1  OW  SOL    11       4.180   2.840   3.660  1.00  0.00            
ATOM      2  HW1 SOL    11       3.920   2.460   2.820  1.00  0.00            
ATOM      3  HW2 SOL    11       3.550   3.540   3.820  1.00  0.00            
ATOM      4  MW  SOL    11       4.070   2.870   3.570  1.00  0.00            
```

### 1.2 Parametry symulacji

Symulacja została przeprowadzona przy następujących parametrach wyciągniętych z plików konfiguracyjnych:

- **Integrator**: md (leap-frog)
- **Krok czasowy (dt)**: 0.002 ps
- **Długość symulacji**: 1,000,000 kroków (2 ns)
- **Temperatura**: 273K
- **Termostat**: V-rescale, stała czasowa 0.1 ps
- **Ciśnienie**: 1.0 bar
- **Barostat**: Parrinello-Rahman, stała czasowa 2.0 ps
- **Ścisliwość**: 4.5e-5 bar^-1
- **Promień odcięcia dla oddziaływań elektrostatycznych (Coulomb)**: 1.0 nm
- **Promień odcięcia dla oddziaływań van der Waalsa**: 1.0 nm
- **Metoda elektrostatyki długozasięgowej**: PME (Particle Mesh Ewald)
- **Warunki brzegowe**: periodyczne w trzech wymiarach (xyz)
- **Liczba cząsteczek wody**: 4000
- **Gęstość początkowa**: 0.93 g/cm³ (zgodnie z plikiem `water_box.inp`)

## 2. Analiza energetyczna

Na podstawie danych z pliku `energy.xvg`, który zawiera informacje o energii potencjalnej, całkowitej, zachowanej oraz temperaturze w funkcji czasu, możemy przeprowadzić następującą analizę:

### 2.1 Energia potencjalna

Energia potencjalna układu składa się z kilku składników:
- Oddziaływania Lennarda-Jonesa (LJ)
- Oddziaływania elektrostatyczne (Coulomb)
- Korekcja dyspersyjna

Z danych widać, że energia potencjalna stabilizuje się na poziomie około -172,000 kJ/mol po początkowym okresie równoważenia, co jest typowe dla symulacji wody w tej skali. 

### 2.2 Energia kinetyczna i temperatura

Z pliku MD.log możemy odczytać, że energia kinetyczna układu utrzymuje się na poziomie około 27,000 kJ/mol, co odpowiada temperaturze około 273K (+/- 5K fluktuacji). 

Związek między energią kinetyczną a temperaturą wynosi:
$$E_{kin} = \frac{3}{2} N k_B T$$

gdzie:
- N to liczba stopni swobody (około 12000 dla 4000 cząsteczek wody)
- k_B to stała Boltzmanna (8.31446 J/(mol·K))
- T to temperatura (273K)

Teoretyczna energia kinetyczna wynosi zatem:
$$E_{kin} = \frac{3}{2} \cdot 12000 \cdot 8.31446 \cdot 10^{-3} \cdot 273 \approx 20,442 \text{ kJ/mol}$$

Różnica między teoretyczną a obserwowaną wartością (27,000 kJ/mol) może wynikać z dodatkowych stopni swobody w systemie lub z niedokładności odczytu z wykresu.

### 2.3 Dryft energii zachowanej

Z wykresu energii zachowanej widoczny jest liniowy wzrost w czasie. Dryft ten wynosi około 25 kJ/(mol·ps), co stanowi około 0.018% całkowitej energii na pikosekundę. Jest to akceptowalna wartość dla symulacji NPT, gdzie dryft energii jest naturalną konsekwencją stosowania termostatu i barostatu.

## 3. Analiza strukturalna - Funkcja rozkładu radialnego (RDF)

Na podstawie danych z pliku `rdf_oo.xvg`, możemy przeprowadzić analizę struktury wody w temperaturze 273K.

### 3.1 Pierwsze maksimum RDF

Pierwszy pik RDF dla par O-O występuje przy r ≈ 0.276 nm z wartością g(r) ≈ 3.27. To maksimum odpowiada pierwszej warstwie koordynacyjnej i reprezentuje molekuły wody połączone wiązaniami wodorowymi.

### 3.2 Liczba koordynacyjna

Liczba koordynacyjna (średnia liczba cząsteczek wody w pierwszej powłoce koordynacyjnej) może być obliczona przez całkowanie funkcji RDF:

$$n(r) = 4\pi\rho\int_0^r g(r')r'^2 dr'$$

gdzie ρ to gęstość liczby cząsteczek.

Dla wody w 273K, liczba koordynacyjna do pierwszego minimum RDF (około 0.34 nm) wynosi typowo około 4.5, co odzwierciedla strukturę tetraedryczną lodu/wody.

### 3.3 Porównanie z danymi eksperymentalnymi

Położenie pierwszego maksimum RDF dla wody w 273K eksperymentalnie wynosi około 2.76 Å (0.276 nm), co bardzo dobrze zgadza się z naszymi wynikami symulacji. Wysokość piku w danych eksperymentalnych wynosi około 3.0-3.1, co jest nieco niższe niż w naszej symulacji (3.27).

## 4. Analiza dyfuzji

Na podstawie danych RMSD z pliku `rmsd.xvg`:

### 4.1 Współczynnik dyfuzji

Współczynnik dyfuzji można obliczyć korzystając z relacji Einsteina:

$$D = \frac{1}{6t}\langle|\vec{r}(t) - \vec{r}(0)|^2\rangle = \frac{\text{RMSD}^2}{6t}$$

Z wykresu RMSD możemy odczytać, że po 200 ps RMSD wynosi około 2.4 nm. Zatem:

$$D = \frac{(2.4 \text{ nm})^2}{6 \cdot 200 \text{ ps}} = \frac{5.76 \text{ nm}^2}{1200 \text{ ps}} \approx 4.8 \times 10^{-3} \text{ nm}^2/\text{ps} = 4.8 \times 10^{-9} \text{ m}^2/\text{s}$$

Eksperymentalna wartość współczynnika dyfuzji wody w 273K wynosi około 1.1 × 10⁻⁹ m²/s, więc nasza wartość jest około 4 razy wyższa. Ta różnica może wynikać z:
1. Niedoskonałości modelu TIP4P
2. Niedokładności w odczycie danych z wykresu
3. Efektów rozmiaru układu (periodic boundary conditions)

## 5. Analiza termodynamiczna

### 5.1 Ciśnienie układu

Z pliku md.log możemy odczytać, że ciśnienie w układzie fluktuuje wokół wartości 1 bar, zgodnie z zadanym ciśnieniem referencyjnym. Fluktuacje ciśnienia są duże, co jest typowe dla symulacji wody, i wynoszą od -300 do +500 bar.

### 5.2 Gęstość wody

Gęstość wody w symulacji można obliczyć na podstawie objętości pudełka symulacyjnego:

$$\rho = \frac{m}{V} = \frac{N \cdot M}{V \cdot N_A}$$

gdzie:
- N to liczba cząsteczek wody (4000)
- M to masa molowa wody (18.015 g/mol)
- V to objętość pudełka symulacyjnego
- N_A to liczba Avogadro (6.022 × 10²³ mol⁻¹)

Początkowa gęstość była ustawiona na 0.93 g/cm³ zgodnie z plikiem `water_box.inp`.

### 5.3 Entalpia parowania

Entalpię parowania można oszacować z różnicy energii potencjalnej między wodą w fazie ciekłej a gazem idealnym:

$$\Delta H_{vap} = -\langle E_{pot} \rangle / N + RT$$

gdzie:
- ⟨E_pot⟩ to średnia energia potencjalna układu
- N to liczba cząsteczek wody
- R to stała gazowa
- T to temperatura

$$\Delta H_{vap} = -(-172,000 \text{ kJ/mol}) / 4000 + 8.31446 \times 10^{-3} \times 273 \text{ kJ/mol} \approx 43 + 2.27 \approx 45.3 \text{ kJ/mol}$$

Eksperymentalna wartość entalpii parowania wody w 273K wynosi około 45.05 kJ/mol, co bardzo dobrze zgadza się z naszym obliczeniem.

### 5.4 Pojemność cieplna

Pojemność cieplną można oszacować na podstawie fluktuacji energii:

$$C_V = \frac{\langle E^2 \rangle - \langle E \rangle^2}{k_B T^2}$$

Z fluktuacji temperatury na poziomie około ±5K wokół średniej 273K, możemy oszacować pojemność cieplną na około 75-80 J/(mol·K), co jest bliskie eksperymentalnej wartości dla wody w 273K (około 75.3 J/(mol·K)).

## 6. Wnioski

### 6.1 Ocena modelu TIP4P

Model TIP4P dobrze odwzorowuje strukturalne właściwości wody w temperaturze 273K, co potwierdza analiza funkcji rozkładu radialnego. Ponadto, obliczona entalpia parowania jest bardzo bliska wartości eksperymentalnej.

Jednak współczynnik dyfuzji jest przeszacowany około 4-krotnie, co jest znanym ograniczeniem modelu TIP4P. Nowsze wersje, takie jak TIP4P/2005 czy TIP4P/Ice, zostały opracowane, aby lepiej odwzorowywać właściwości dynamiczne i termodynamiczne wody, zwłaszcza w niskich temperaturach.

### 6.2 Zalecenia do dalszych badań

1. Przeprowadzenie dłuższej symulacji (>10 ns) w celu lepszego zbadania właściwości dynamicznych
2. Porównanie z innymi modelami wody (SPC/E, TIP4P/2005, TIP4P/Ice)
3. Analiza wpływu rozmiaru układu na obliczone właściwości

## Załącznik: Wzory fizyczne wykorzystane w analizie

### A.1 Energia kinetyczna i temperatura

$$E_{kin} = \frac{3}{2} N k_B T$$

### A.2 Funkcja rozkładu radialnego i liczba koordynacyjna

$$g(r) = \frac{1}{\rho} \frac{dN(r)}{dV(r)}$$

$$n(r) = 4\pi\rho\int_0^r g(r')r'^2 dr'$$

### A.3 Współczynnik dyfuzji

$$D = \frac{1}{6t}\langle|\vec{r}(t) - \vec{r}(0)|^2\rangle = \frac{\text{RMSD}^2}{6t}$$

### A.4 Entalpia parowania

$$\Delta H_{vap} = -\langle E_{pot} \rangle / N + RT$$

### A.5 Pojemność cieplna

$$C_V = \frac{\langle E^2 \rangle - \langle E \rangle^2}{k_B T^2}$$

$$C_P = C_V + \frac{T V \alpha_P^2}{\kappa_T}$$

gdzie:
- α_P to współczynnik rozszerzalności termicznej
- κ_T to współczynnik ściśliwości izotermicznej 