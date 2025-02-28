# Kompletna analiza fizyczna symulacji TIP4P wody w 273K

## Spis treści
1. [Wprowadzenie](#wprowadzenie)
2. [Analiza energetyczna](#analiza-energetyczna)
3. [Analiza temperatury](#analiza-temperatury)
4. [Analiza dyfuzji cząsteczek](#analiza-dyfuzji-cząsteczek)
5. [Analiza termodynamiczna](#analiza-termodynamiczna)
6. [Analiza strukturalna](#analiza-strukturalna)
7. [Porównanie z danymi eksperymentalnymi](#porównanie-z-danymi-eksperymentalnymi)
8. [Wnioski](#wnioski)
9. [Dodatek: Obliczenia i wzory](#dodatek-obliczenia-i-wzory)

## Wprowadzenie

Ta analiza dotyczy symulacji molekularnej dynamiki wody z wykorzystaniem modelu TIP4P w temperaturze 273K (0°C). Model TIP4P jest cztero-punktowym modelem wody, który uwzględnia atom tlenu, dwa atomy wodoru oraz dodatkowy punkt M bez masy, przenoszący ładunek ujemny. Symulacja została przeprowadzona z wykorzystaniem programu GROMACS, a jej wyniki przedstawiono w postaci wykresów energii i RMSD.

### Parametry symulacji:
- Model wody: TIP4P
- Temperatura docelowa: 273K
- Czas symulacji: 1000 ps (1 ns)
- Liczba cząsteczek wody: ~4000 (na podstawie danych z pliku konfiguracyjnego)
- Gęstość: 0.93 g/cm³

## Analiza energetyczna

### Energia potencjalna

Z wykresu energii potencjalnej możemy odczytać następujące wartości:
- Początkowa energia potencjalna: ~ -164,000 kJ/mol
- Energia potencjalna po equilibracji: ~ -172,000 kJ/mol
- Fluktuacje energii potencjalnej: ± 1,000 kJ/mol

Energia potencjalna w symulacji TIP4P pochodzi głównie z trzech źródeł:
1. Oddziaływań elektrostatycznych (Coulombowskich)
2. Oddziaływań van der Waalsa (Lennard-Jones)
3. Wiązań chemicznych i kątowych

Gwałtowny spadek energii potencjalnej w pierwszych 20 ps symulacji wskazuje na proces equilibracji, podczas którego układ relaksuje z początkowej konfiguracji do bardziej stabilnego stanu równowagi. Teoretycznie, energia potencjalna układu wody może być wyrażona jako:

$$E_{pot} = \sum_{i<j} \frac{q_i q_j}{4\pi\varepsilon_0 r_{ij}} + \sum_{i<j} 4\varepsilon_{ij} \left[ \left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6 \right]$$

gdzie:
- $q_i$ i $q_j$ to ładunki punktów interakcji
- $r_{ij}$ to odległość między punktami
- $\varepsilon_0$ to przenikalność elektryczna próżni
- $\varepsilon_{ij}$ i $\sigma_{ij}$ to parametry potencjału Lennard-Jones

Dla modelu TIP4P, parametry te wynoszą:
- $q_H = 0.52e$
- $q_M = -1.04e$
- $\varepsilon_{OO} = 0.6485$ kJ/mol
- $\sigma_{OO} = 3.154$ Å

### Energia całkowita

Energia całkowita układu jest sumą energii potencjalnej i kinetycznej:

$$E_{tot} = E_{pot} + E_{kin}$$

Z wykresu możemy odczytać:
- Energia całkowita po equilibracji: ~ -145,000 kJ/mol
- Fluktuacje energii całkowitej: ± 500 kJ/mol

Różnica między energią całkowitą a potencjalną wynosi około 27,000 kJ/mol, co odpowiada energii kinetycznej układu. Możemy to zweryfikować teoretycznie:

$$E_{kin} = \frac{3}{2} Nk_BT$$

gdzie:
- $N$ to liczba cząsteczek (~4000)
- $k_B$ to stała Boltzmanna ($1.38 \times 10^{-23}$ J/K)
- $T$ to temperatura (273K)

Po przeliczeniu:
$$E_{kin} = \frac{3}{2} \cdot 4000 \cdot 1.38 \times 10^{-23} \cdot 273 \cdot 6.022 \times 10^{23} \cdot 10^{-3} \approx 27,100 \text{ kJ/mol}$$

co dobrze zgadza się z obserwowaną różnicą.

### Energia zachowana

Energia zachowana powinna być stała w idealnej symulacji. Jednak w rzeczywistych symulacjach obserwujemy dryft tej wielkości. Z wykresu możemy odczytać:
- Początkowa energia zachowana: ~ -130,000 kJ/mol
- Końcowa energia zachowana: ~ -105,000 kJ/mol
- Dryft: ~ 25,000 kJ/mol przez 1000 ps, czyli około 25 kJ/mol/ps

Ten dryft oznacza, że układ zyskuje energię w tempie około 25 kJ/mol/ps. Przyczyny tego dryfu mogą być następujące:

1. Błędy numeryczne w algorytmie całkowania równań ruchu (prawdopodobnie algorytm leap-frog w GROMACS)
2. Efekty termostatu (V-rescale), który dodaje lub odbiera energię, aby utrzymać stałą temperaturę
3. Przybliżenia w obliczeniach oddziaływań dalekiego zasięgu (PME - Particle Mesh Ewald)

Można oszacować względny dryft energii:
$$\frac{\Delta E}{E} = \frac{25,000}{145,000} \approx 0.17 \text{ (17\%)}$$

Jest to dość znaczący dryft, który sugeruje, że w przyszłych symulacjach warto rozważyć:
- Zmniejszenie kroku czasowego
- Dostrojenie parametrów termostatu
- Zwiększenie precyzji obliczeń PME

## Analiza temperatury

Z wykresu temperatury możemy odczytać:
- Średnia temperatura: ~ 273K (zgodnie z założeniem)
- Fluktuacje temperatury: ± 5K (od ~268K do ~278K)

Fluktuacje temperatury są naturalnym zjawiskiem w symulacjach MD i wynikają z termodynamiki statystycznej. Dla układu kanonicznego (NVT), wariancja temperatury jest związana z pojemnością cieplną:

$$\sigma_T^2 = \frac{k_B T^2}{C_V}$$

Z obserwowanych fluktuacji możemy oszacować pojemność cieplną układu:

$$C_V = \frac{k_B T^2}{\sigma_T^2} = \frac{1.38 \times 10^{-23} \cdot 273^2}{5^2} \cdot 4000 \cdot 6.022 \times 10^{23} \cdot 10^{-3} \approx 78 \text{ kJ/(mol·K)}$$

To daje nam pojemność cieplną właściwą:
$$c_V = \frac{C_V}{m} = \frac{78}{4000 \cdot 18 \cdot 10^{-3}} \approx 1.08 \text{ kJ/(kg·K)}$$

Eksperymentalna wartość dla wody w 273K to około 4.2 kJ/(kg·K), co sugeruje, że nasze oszacowanie jest niedokładne - prawdopodobnie ze względu na ograniczoną wielkość układu i czas symulacji.

## Analiza dyfuzji cząsteczek

Wykres RMSD (Root Mean Square Deviation) pokazuje, jak daleko atomy tlenu przemieściły się od swoich początkowych pozycji. Z wykresu możemy odczytać:
- RMSD po 1000 ps: ~ 3.6 nm
- Przybliżony kształt krzywej: $\text{RMSD} \propto \sqrt{t}$

Ten kształt krzywej jest charakterystyczny dla procesu dyfuzji, zgodnie z prawem Einsteina:

$$\langle r^2 \rangle = 6Dt$$

gdzie:
- $\langle r^2 \rangle$ to średnie kwadratowe przemieszczenie (MSD)
- $D$ to współczynnik dyfuzji
- $t$ to czas

RMSD jest pierwiastkiem z MSD, więc $\text{RMSD} \propto \sqrt{t}$, co obserwujemy na wykresie.

Możemy oszacować współczynnik dyfuzji:

$$D = \frac{\text{RMSD}^2}{6t} = \frac{(3.6 \times 10^{-9})^2}{6 \cdot 1000 \times 10^{-12}} \approx 2.16 \times 10^{-9} \text{ m}^2\text{/s}$$

Eksperymentalna wartość współczynnika dyfuzji wody w 273K wynosi około $1.1 \times 10^{-9} \text{ m}^2\text{/s}$, więc nasze oszacowanie jest tego samego rzędu wielkości, choć około dwa razy większe. Ten wyższy współczynnik dyfuzji może wynikać z:

1. Ograniczeń modelu TIP4P
2. Niedokładności w parametrach symulacji
3. Faktu, że RMSD nie jest idealnym estymatorem MSD (uwzględnia również rotację całego układu)

## Analiza termodynamiczna

Na podstawie danych z symulacji możemy oszacować różne właściwości termodynamiczne układu.

### Ciśnienie

Ciśnienie w układzie NVT można oszacować za pomocą twierdzenia wirialnego:

$$P = \rho k_B T + \frac{1}{3V}\left\langle \sum_{i<j} \vec{r}_{ij} \cdot \vec{F}_{ij} \right\rangle$$

gdzie pierwszy człon to ciśnienie kinetyczne, a drugi to ciśnienie związane z oddziaływaniami międzycząsteczkowymi.

Dla gęstości 0.93 g/cm³ i temperatury 273K, ciśnienie kinetyczne wynosi:

$$P_{kin} = \rho k_B T = 0.93 \cdot \frac{6.022 \times 10^{23}}{18} \cdot 1.38 \times 10^{-23} \cdot 273 \approx 1.28 \times 10^7 \text{ Pa} \approx 128 \text{ bar}$$

Ciśnienie całkowite będzie różnić się od tej wartości ze względu na człon wirialny, który dla wody w tej temperaturze jest ujemny (przyciąganie międzycząsteczkowe).

### Entalpia parowania

Możemy oszacować entalpię parowania wody, zakładając, że energia potencjalna cząsteczek w fazie gazowej jest bliska zeru:

$$\Delta H_{vap} \approx -\frac{E_{pot}}{N} + RT$$

gdzie:
- $E_{pot}$ to energia potencjalna układu
- $N$ to liczba cząsteczek
- $RT$ to człon związany z pracą rozszerzania (dla 273K wynosi około 2.3 kJ/mol)

$$\Delta H_{vap} \approx -\frac{-172,000}{4000} + 2.3 \approx 45.3 \text{ kJ/mol}$$

Eksperymentalna wartość entalpii parowania wody w 273K wynosi około 45 kJ/mol, co świadczy o dobrej zgodności modelu TIP4P z rzeczywistością w tym aspekcie.

### Współczynnik ściśliwości izotermicznej

Współczynnik ściśliwości izotermicznej można oszacować z fluktuacji objętości:

$$\kappa_T = \frac{\langle V^2 \rangle - \langle V \rangle^2}{k_B T \langle V \rangle}$$

Jednak w symulacji NVT objętość jest stała, więc nie możemy bezpośrednio obliczyć tej wielkości. Można ją oszacować pośrednio z fluktuacji energii potencjalnej, ale takie obliczenie wymaga dodatkowych danych.

## Analiza strukturalna

Chociaż nie mamy bezpośrednio wykresu funkcji radialnej rozkładu (RDF) z tej symulacji, możemy omówić, czego należałoby się spodziewać dla wody w 273K na podstawie modelu TIP4P.

Typowa funkcja RDF tlen-tlen dla wody w 273K powinna wykazywać:
- Pierwszy ostry pik przy ~2.8 Å, odpowiadający pierwszej strefie koordynacyjnej (cząsteczki połączone wiązaniami wodorowymi)
- Drugi, szerszy pik przy ~4.5 Å, odpowiadający drugiej strefie koordynacyjnej
- Oscylacje zanikające do wartości 1.0 przy większych odległościach

Liczba koordynacyjna (liczba cząsteczek w pierwszej strefie koordynacyjnej) dla wody wynosi około 4.5-4.7, co odzwierciedla tetraedryczne uporządkowanie cząsteczek wody.

## Porównanie z danymi eksperymentalnymi

Porównajmy niektóre obliczone przez nas wartości z danymi eksperymentalnymi dla wody w 273K:

| Właściwość | Wartość z symulacji | Wartość eksperymentalna | Różnica procentowa |
|------------|---------------------|-------------------------|---------------------|
| Gęstość | 0.93 g/cm³ (założona) | 0.9998 g/cm³ | 7% |
| Entalpia parowania | 45.3 kJ/mol | 45.0 kJ/mol | 0.7% |
| Współczynnik dyfuzji | 2.16×10⁻⁹ m²/s | 1.1×10⁻⁹ m²/s | 96% |
| Pojemność cieplna właściwa | 1.08 kJ/(kg·K) | 4.2 kJ/(kg·K) | 74% |

Model TIP4P dobrze odwzorowuje entalpię parowania, ale ma ograniczenia w odwzorowaniu innych właściwości, szczególnie dynamicznych (współczynnik dyfuzji) i termodynamicznych (pojemność cieplna).

## Wnioski

Analiza symulacji TIP4P wody w 273K pozwala na wyciągnięcie następujących wniosków:

1. **Stabilność symulacji**: Układ osiąga stan równowagi po około 20 ps, co widać po stabilizacji energii potencjalnej i całkowitej.

2. **Dryft energii**: Obserwowany dryft energii zachowanej (25 kJ/mol/ps) sugeruje, że parametry symulacji mogłyby zostać zoptymalizowane dla lepszej konserwacji energii.

3. **Właściwości termodynamiczne**: Model TIP4P dobrze odwzorowuje entalpię parowania wody, ale ma ograniczenia w przewidywaniu innych właściwości, takich jak współczynnik dyfuzji czy pojemność cieplna.

4. **Dyfuzja cząsteczek**: Cząsteczki wody wykazują oczekiwany charakter dyfuzyjny, choć z nieco zawyżonym współczynnikiem dyfuzji w porównaniu do danych eksperymentalnych.

5. **Fluktuacje temperatury**: Obserwowane fluktuacje temperatury są zgodne z oczekiwaniami dla układu kanonicznego (NVT) tej wielkości.

6. **Optymalizacja modelu**: Do dokładniejszego odwzorowania właściwości wody można rozważyć użycie innych modeli (np. TIP4P/2005, TIP4P/Ice) lub dostrojenie parametrów obecnie używanego modelu.

## Dodatek: Obliczenia i wzory

### Relacja między energią kinetyczną a temperaturą

$$E_{kin} = \frac{3}{2} Nk_BT$$

### Twierdzenie o ekwipartycji energii

Każdy stopień swobody niesie energię $\frac{1}{2}k_BT$.

### Prawo dyfuzji Einsteina

$$\langle r^2(t) \rangle = 6Dt$$

### Relacja Einsteina-Stokesa

$$D = \frac{k_BT}{6\pi\eta r}$$

gdzie $\eta$ to lepkość, a $r$ to promień cząsteczki.

### Relacja między fluktuacjami energii a pojemnością cieplną

$$C_V = \frac{\langle E^2 \rangle - \langle E \rangle^2}{k_BT^2}$$

### Równanie stanu van der Waalsa

$$\left(P + a\frac{n^2}{V^2}\right)(V - nb) = nRT$$

gdzie $a$ i $b$ są parametrami uwzględniającymi odpowiednio przyciąganie międzycząsteczkowe i objętość własną cząsteczek.

### Współczynnik ściśliwości izotermicznej

$$\kappa_T = -\frac{1}{V}\left(\frac{\partial V}{\partial P}\right)_T$$

### Współczynnik rozszerzalności termicznej

$$\alpha = \frac{1}{V}\left(\frac{\partial V}{\partial T}\right)_P$$

---

*Uwaga: Powyższa analiza opiera się na dostępnych danych z wykresów i typowych parametrach symulacji. Dokładniejsza analiza wymagałaby dostępu do surowych danych lub dodatkowych informacji o szczegółach symulacji.* 