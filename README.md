# Spectrum Analysis Tool

Narzędzie do analizy spektralnej danych autokorelacyjnych.

## Funkcje:
- Symetryzacja funkcji autokorelacyjnej.
- Usuwanie składowej stałej (DC offset).
- Obliczanie FFT i widma mocy za pomocą `pyFFTW`.
- Wybór różnych typów transformacji (FFT, rFFT, DCT).
- Generowanie dynamicznego raportu HTML z wynikami.
- Obsługa błędów przy wczytywaniu plików i analizie.

## Instalacja

1. Sklonuj repozytorium:
    ```bash
    git clone https://github.com/TWOJ_USERNAME/spectrum_analysis_tool.git
    cd spectrum_analysis_tool
    ```

2. Zainstaluj wymagane pakiety:
    ```bash
    pip install -r requirements.txt
    ```

## Konfiguracja

Skonfiguruj narzędzie za pomocą pliku `config.yaml` (przykład poniżej).

```yaml
data_dir: 'data/'
output_dir: 'output/'
sampling_frequency: 1000
num_functions: 10
samples_per_function: 1024
log_file: 'spectrum_analysis.log'
log_level: 'INFO'
fft_type: 'rfft'
report_file: 'raport.html'




Spectrum Analysis
Autorzy
Wojciech Żurawski – GitHub | watwzwp@gmail.com
Eugeniusz Pazderski
Spis treści
Opis
Wymagania
Instalacja
Struktura projektu
Przygotowanie danych
Format danych wejściowych
Plik wykluczeń
Użycie
Przykładowe wyniki
Raport
Licencja
Dodatkowe Uwagi
Opis
Spectrum Analysis to program napisany w Pythonie, który analizuje funkcje autokorelacyjne i oblicza widma mocy za pomocą szybkiej transformaty Fouriera (FFT). Program jest w stanie przetwarzać wiele plików danych jednocześnie, wykluczając określone pliki na podstawie listy wykluczeń. Wyniki analizy są przedstawiane w formie wykresów statycznych i interaktywnych oraz zawarte w szczegółowym raporcie HTML.

Wymagania
Python 3.6 lub nowszy
Pakiety wymienione w requirements.txt:
numpy
scipy
matplotlib
plotly
pyyaml

1. Instalacja
Sklonuj repozytorium:
git clone https://github.com/wzurawski015/spectrum_analysis.git
cd spectrum_analysis

2. Stwórz środowisko wirtualne (opcjonalnie):
python3 -m venv venv
source venv/bin/activate  # Na Windows: venv\Scripts\activate

3. Zainstaluj zależności:
pip install -r requirements.txt

Struktura projektu
spectrum_analysis/
├── spectrum_analysis.py      # Główny skrypt programu
├── requirements.txt          # Plik z zależnościami
├── README.md                 # Dokumentacja projektu
├── LICENSE                   # Informacje o licencji
├── .gitignore                # Plik konfiguracji Git
├── config.yaml               # Plik konfiguracyjny
├── data/                     # Folder z danymi wejściowymi
│   ├── input1.dat            # Przykładowy plik danych
│   ├── input2.dat            # Kolejny plik danych
│   ├── sample3.dat           # Inny plik danych
│   └── exclude               # Plik z nazwami plików do wykluczenia (opcjonalny)
└── output/                   # Folder z wynikami (generowany automatycznie)
    ├── input1_pcal1_autocorr.png          # Wykres funkcji autokorelacyjnej 1
    ├── input1_pcal1_fft.png                # Wykres widma mocy 1
    ├── input1_pcal1_interactive.html       # Interaktywny wykres widma mocy 1
    ├── input1_pcal1_intermediate.json      # Pośrednie wyniki obliczeń 1
    ├── raport.html                         # Raport HTML z wynikami analizy
    └── ...                                 # Inne wyniki analizy

Przygotowanie danych
Format danych wejściowych
Pliki z danymi wejściowymi powinny zawierać cztery funkcje autokorelacyjne, każda reprezentowana przez pary:
Nr_kanalu  wartosc_funkcji

Nr_kanalu: Numer kanału (lag/opóźnienie), liczba całkowita.
wartosc_funkcji: Wartość funkcji autokorelacyjnej dla danego numeru kanału, liczba zmiennoprzecinkowa.
Wymagania dotyczące danych:

Numery kanałów powinny być od 0 do 16387.
Dane dla każdej funkcji powinny być umieszczone kolejno, bez przerw między nimi.
Każda funkcja autokorelacyjna powinna zawierać dokładnie 4097 próbek dla funkcji 1-3 i 4096 próbek dla funkcji 4.
Przykład fragmentu pliku input1.dat:
0       124413044
1       95279359
2       64397895
...
16387   62230859

Plik wykluczeń
Aby wykluczyć określone pliki z analizy, utwórz plik exclude w katalogu data i dodaj do niego nazwy plików, które mają być pominięte. Każda nazwa pliku powinna być na osobnej linii.

Przykład zawartości pliku data/exclude:
input2.dat
sample3.dat
old_data.dat

Uwagi:

Nazwy plików muszą dokładnie odpowiadać nazwom plików w katalogu data.
Nie uwzględniaj ścieżek, tylko same nazwy plików.
Użycie
1. Przygotuj dane wejściowe:
Umieść wszystkie pliki z danymi wejściowymi w folderze data/.
Jeśli chcesz wykluczyć niektóre pliki, utwórz plik data/exclude z listą nazw plików do pominięcia.
2. Uruchom program:
python spectrum_analysis.py --config config.yaml

3. Sprawdź wyniki:

Wyniki analizy zostaną zapisane w folderze output/, w tym:

Pośrednie pliki z obliczeniami (*_intermediate.json): Zawierają przetworzone dane autokorelacyjne, częstotliwości FFT oraz widmo mocy.
Wykresy funkcji autokorelacyjnych (*_autocorr.png): Przedstawiają graficzną interpretację funkcji autokorelacyjnych.
Wykresy widma mocy (*_fft.png): Przedstawiają graficzną interpretację widm mocy.
Interaktywne wykresy widma mocy (*_interactive.html): Umożliwiają interaktywną eksplorację widma mocy w przeglądarce internetowej.
Raport HTML (raport.html): Zawiera wszystkie powyższe wyniki w jednym, przejrzystym pliku HTML

Przykładowe wyniki
Po uruchomieniu programu w folderze output/ znajdziesz:

Pliki z funkcjami autokorelacyjnymi (input1_pcal1_intermediate.json, ...):

Zawierają przetworzone i symetryzowane funkcje autokorelacyjne, częstotliwości FFT oraz widma mocy w formacie JSON.
Pliki z widmami mocy (input1_pcal1_fft.png, ...):

Zawierają wykresy widm mocy w skali decybelowej (dB).
Wykresy funkcji autokorelacyjnych (input1_pcal1_autocorr.png, ...):

Przedstawiają graficzną interpretację funkcji autokorelacyjnych.
Interaktywne wykresy widma mocy (input1_pcal1_interactive.html, ...):

Umożliwiają interaktywną eksplorację widma mocy w przeglądarce internetowej.
Raport HTML (raport.html):

Zawiera wszystkie powyższe wyniki w jednym, przejrzystym pliku HTML.

Raport
Raport HTML (raport.html) zawiera:

Nagłówek z tytułem i datą generacji raportu.
Sekcje dla każdego przetworzonego pliku danych:
Nazwa pliku danych.
Wykres funkcji autokorelacyjnej.
Wykres widma mocy.
Link do interaktywnego wykresu widma mocy.
Przykład sekcji dla jednego pliku danych:
<div class="file-section">
    <h2>Plik danych: input1.dat</h2>
    <h3>Funkcja autokorelacyjna 1</h3>
    <p><strong>Funkcja autokorelacyjna 1:</strong></p>
    <img src="input1_pcal1_autocorr.png" alt="Autokorelacja 1" width="800">
    <p><strong>Widmo mocy 1:</strong></p>
    <img src="input1_pcal1_fft.png" alt="Widmo mocy 1" width="800">
    <p><strong>Interaktywny wykres widma mocy 1:</strong> <a href="input1_pcal1_interactive.html" target="_blank">Otwórz</a></p>
    <hr>
</div>

Licencja
Ten projekt jest objęty licencją MIT License – więcej informacji w pliku LICENSE.

Autorzy
Wojciech Żurawski – GitHub | watwzwp@gmail.com
Eugeniusz Pazderski
Dodatkowe Uwagi
Częstotliwość próbkowania (FS)
Upewnij się, że wartość FS w pliku config.yaml odpowiada częstotliwości próbkowania Twoich danych. Domyślnie jest ustawiona na 1000 Hz. Jeśli Twoje dane mają inną częstotliwość próbkowania, zmień tę wartość na odpowiednią.
FS: 1000  # Częstotliwość próbkowania w Hz

Skrypt jest skonfigurowany do analizowania wszystkich plików w katalogu data/ niezależnie od ich rozszerzenia. Jeśli chcesz ograniczyć analizę do określonych rozszerzeń (np. .dat, .txt), zmodyfikuj linię w funkcji main():
data_files = glob.glob(os.path.join(DATA_DIR, '*'))

Na przykład, aby analizować tylko pliki .dat:
data_files = glob.glob(os.path.join(DATA_DIR, '*.dat'))

Obsługa błędów
Skrypt informuje użytkownika o wszelkich problemach z danymi, takich jak nieprawidłowa liczba wierszy czy błędy podczas zapisywania plików. Ułatwia to diagnozowanie i rozwiązywanie problemów.

Dodawanie dodatkowych funkcjonalności
Jeśli chcesz dodać więcej szczegółów do raportu, takich jak statystyki opisowe, histogramy czy inne rodzaje wykresów, możesz rozszerzyć funkcję generate_report() oraz dodać odpowiednie funkcje do przetwarzania i wizualizacji danych.

Podsumowanie
Spectrum Analysis to wszechstronny narzędzie do analizy spektralnej danych autokorelacyjnych, które umożliwia przetwarzanie wielu plików jednocześnie, generowanie różnorodnych wykresów oraz tworzenie szczegółowego raportu HTML. Dzięki przetwarzaniu równoległemu, program jest wydajny nawet przy dużej liczbie plików danych. Dzięki plikowi konfiguracyjnemu config.yaml, użytkownicy mogą łatwo dostosować parametry analizy do swoich potrzeb.

Dziękujemy za skorzystanie z Spectrum Analysis! Jeśli masz pytania lub sugestie dotyczące projektu, zapraszamy do kontaktu poprzez zgłoszenie problemu (issue) na GitHubie lub bezpośredni kontakt z autorami.
