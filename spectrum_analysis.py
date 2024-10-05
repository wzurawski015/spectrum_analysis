#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spectrum_analysis.py

Skrypt do analizy spektralnej danych autokorelacyjnych.
Wykonuje następujące kroki:
1. Wczytuje dane z plików tekstowych.
2. Ignoruje pliki na liście wykluczeń.
3. Przetwarza każdą funkcję autokorelacyjną:
   a. Symetryzuje funkcję autokorelacyjną.
   b. Usuwa składową stałą (DC offset).
   c. Oblicza FFT i widmo mocy.
   d. Zapisuje pośrednie wyniki do plików JSON.
   e. Generuje wykresy statyczne i interaktywne.
4. Tworzy raport HTML zawierający wszystkie wyniki analizy.
5. Wykorzystuje przetwarzanie równoległe dla zwiększenia wydajności.

Autor: [Twoje Imię]
Data: 2024-04-27
"""

import os
import glob
import argparse
import logging
import json
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import yaml
from scipy.fft import fft, rfft, fftfreq, rfftfreq
from datetime import datetime
from multiprocessing import Pool, cpu_count

def load_config(config_path):
    """
    Wczytuje konfigurację z pliku YAML.

    Args:
        config_path (str): Ścieżka do pliku konfiguracyjnego YAML.

    Returns:
        dict: Słownik zawierający konfigurację.
    """
    with open(config_path, 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)
    return config

def setup_logging(log_file, log_level):
    """
    Konfiguruje logowanie do pliku i konsoli.

    Args:
        log_file (str): Ścieżka do pliku logu.
        log_level (str): Poziom logowania (np. INFO, DEBUG).
    """
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        numeric_level = logging.INFO
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, encoding='utf-8'),
            logging.StreamHandler()
        ]
    )

def parse_arguments():
    """
    Parsuje argumenty wiersza poleceń.

    Returns:
        argparse.Namespace: Obiekt zawierający argumenty.
    """
    parser = argparse.ArgumentParser(description='Spectrum Analysis EP')
    parser.add_argument('--config', type=str, default='config.yaml',
                        help='Ścieżka do pliku konfiguracyjnego')
    return parser.parse_args()

def load_exclude_list(exclude_path):
    """
    Wczytuje listę plików do wykluczenia z pliku exclude.
    Jeśli plik nie istnieje, zwraca pustą listę.

    Args:
        exclude_path (str): Ścieżka do pliku z listą wykluczeń.

    Returns:
        list: Lista nazw plików do wykluczenia.
    """
    exclude_files = []
    if os.path.exists(exclude_path):
        with open(exclude_path, 'r', encoding='utf-8') as f:
            exclude_files = [line.strip() for line in f if line.strip()]
        logging.info(f"Załadowano {len(exclude_files)} plików do wykluczenia.")
    else:
        logging.info("Plik 'exclude' nie istnieje. Analizowane będą wszystkie pliki w katalogu 'data'.")
    return exclude_files

def get_data_files(data_dir, exclude_files):
    """
    Zwraca listę plików do analizy, pomijając te na liście wykluczeń oraz sam plik 'exclude'.

    Args:
        data_dir (str): Katalog z danymi wejściowymi.
        exclude_files (list): Lista plików do wykluczenia.

    Returns:
        list: Lista ścieżek do plików danych do analizy.
    """
    all_files = glob.glob(os.path.join(data_dir, '*'))
    data_files = [
        f for f in all_files
        if os.path.isfile(f)
        and os.path.basename(f) not in exclude_files
        and os.path.basename(f) != 'exclude'
    ]
    logging.info(f"Znaleziono {len(data_files)} plików do analizy.")
    return data_files

def read_autocorrelation(file_path, num_functions, samples_per_function):
    """
    Wczytuje funkcje autokorelacyjne z pliku danych.

    Zakłada, że dane są podzielone na `num_functions` funkcji,
    każda zawierająca `samples_per_function` próbek.

    Args:
        file_path (str): Ścieżka do pliku danych.
        num_functions (int): Liczba funkcji autokorelacyjnych w pliku.
        samples_per_function (int): Liczba próbek na funkcję autokorelacyjną.

    Returns:
        list: Lista funkcji autokorelacyjnych jako numpy arrays.
    """
    try:
        data = np.loadtxt(file_path, delimiter=None)  # Zakładamy, że dane są oddzielone białymi znakami
        if data.ndim == 1:
            data = data.reshape(-1, 2)  # Zakładaj, że dane mają dwie kolumny
        if data.shape[1] < 2:
            raise ValueError("Plik danych powinien mieć co najmniej dwie kolumny.")
        channels = data[:, 0].astype(int)
        values = data[:, 1].astype(float)
        autocorr_functions = []
        total_samples = num_functions * samples_per_function
        if len(values) < total_samples:
            raise ValueError(f"Oczekiwano co najmniej {total_samples} próbek, ale znaleziono {len(values)}.")
        for i in range(num_functions):
            start = i * samples_per_function
            end = start + samples_per_function
            autocorr = values[start:end]
            autocorr_functions.append(autocorr)
        logging.debug(f"Plik {file_path} zawiera {len(autocorr_functions)} funkcji autokorelacyjnych.")
        return autocorr_functions
    except Exception as e:
        logging.error(f"Błąd podczas wczytywania pliku {file_path}: {e}")
        return []

def symmetrize(autocorr):
    """
    Symetryzuje funkcję autokorelacyjną, aby uzyskać rzeczywiste widmo mocy.

    Args:
        autocorr (numpy.ndarray): Funkcja autokorelacyjna.

    Returns:
        numpy.ndarray: Symetryzowana funkcja autokorelacyjna.
    """
    return np.concatenate((autocorr, autocorr[::-1]))

def remove_dc_offset(autocorr):
    """
    Usuwa składową stałą (DC offset) z funkcji autokorelacyjnej.

    Args:
        autocorr (numpy.ndarray): Funkcja autokorelacyjna.

    Returns:
        numpy.ndarray: Funkcja autokorelacyjna bez składowej stałej.
    """
    return autocorr - np.mean(autocorr)

def compute_fft(autocorr, fft_type, fs):
    """
    Wykonuje szybką transformatę Fouriera (FFT) na funkcji autokorelacyjnej.

    Args:
        autocorr (numpy.ndarray): Funkcja autokorelacyjna.
        fft_type (str): Typ FFT do użycia ('fft' lub 'rfft').
        fs (int): Częstotliwość próbkowania w Hz.

    Returns:
        tuple: FFT wynik (numpy.ndarray), częstotliwości (numpy.ndarray).
    """
    if fft_type == 'rfft':
        fft_result = rfft(autocorr)
        freqs = rfftfreq(len(autocorr), d=1/fs)
    else:
        fft_result = fft(autocorr)
        freqs = fftfreq(len(autocorr), d=1/fs)
    return fft_result, freqs

def compute_power_spectrum(fft_result):
    """
    Oblicza widmo mocy i przekształca je na skalę decybelową (dB).

    Args:
        fft_result (numpy.ndarray): Wynik FFT.

    Returns:
        numpy.ndarray: Widmo mocy w skali dB.
    """
    power_spectrum = 10 * np.log10(np.abs(fft_result) + 1e-12)  # Dodanie 1e-12 dla stabilności logarytmu
    return power_spectrum

def save_intermediate_results(file_name, idx, cleaned_autocorr, freqs, power_spectrum, output_dir):
    """
    Zapisuje pośrednie wyniki obliczeń do plików JSON.

    Plik JSON zawiera:
    - cleaned_autocorr: Oczyszczona funkcja autokorelacyjna po symetryzacji i usunięciu DC offset.
    - freqs: Częstotliwości FFT.
    - power_spectrum: Widmo mocy w skali dB.

    Args:
        file_name (str): Nazwa pliku danych.
        idx (int): Indeks funkcji autokorelacyjnej.
        cleaned_autocorr (numpy.ndarray): Oczyszczona funkcja autokorelacyjna.
        freqs (numpy.ndarray): Częstotliwości FFT.
        power_spectrum (numpy.ndarray): Widmo mocy.
        output_dir (str): Katalog z wynikami.
    """
    intermediate_data = {
        'cleaned_autocorr': cleaned_autocorr.tolist(),
        'freqs': freqs.tolist(),
        'power_spectrum': power_spectrum.tolist()
    }
    output_prefix = os.path.join(
        output_dir, f"{os.path.splitext(file_name)[0]}_pcal{idx}"
    )
    json_path = f"{output_prefix}_intermediate.json"
    try:
        with open(json_path, 'w', encoding='utf-8') as json_file:
            json.dump(intermediate_data, json_file, ensure_ascii=False, indent=4)
        logging.debug(f"Pośrednie wyniki zapisane jako {json_path}")
    except Exception as e:
        logging.error(f"Błąd podczas zapisywania pośrednich wyników do {json_path}: {e}")

def generate_plots(args):
    """
    Generuje wykresy funkcji autokorelacyjnych i widm mocy.
    Tworzy zarówno wykresy statyczne (.png), jak i interaktywne (.html).
    Zapisuje również pośrednie wyniki do plików JSON.

    Args:
        args (tuple): Zawiera wszystkie niezbędne argumenty do przetworzenia pliku.

    Returns:
        dict or None: Słownik z informacjami o plikach wykresów, lub None w przypadku błędu.
    """
    file_path, output_dir, fs, num_functions, samples_per_function, fft_type = args
    file_name = os.path.basename(file_path)
    logging.info(f"Przetwarzanie pliku: {file_name}")
    autocorr_functions = read_autocorrelation(file_path, num_functions, samples_per_function)
    if not autocorr_functions:
        logging.warning(f"Pominięto plik {file_name} z powodu błędu podczas wczytywania.")
        return None

    plots_info = []
    for idx, autocorr in enumerate(autocorr_functions, start=1):
        try:
            # Symetryzacja funkcji autokorelacyjnej
            sym_autocorr = symmetrize(autocorr)
            # Usunięcie składowej stałej (DC offset)
            cleaned_autocorr = remove_dc_offset(sym_autocorr)
            # Obliczenie FFT
            fft_result, freqs = compute_fft(cleaned_autocorr, fft_type, fs)
            # Obliczenie widma mocy
            power_spectrum = compute_power_spectrum(fft_result)
            # Zapisanie pośrednich wyników
            save_intermediate_results(file_name, idx, cleaned_autocorr, freqs, power_spectrum, output_dir)

            # Generowanie wykresu funkcji autokorelacyjnej
            plt.figure(figsize=(10, 6))
            plt.plot(cleaned_autocorr, label='Autokorelacja')
            plt.title(f'Funkcja Autokorelacyjna {idx} - {file_name}')
            plt.xlabel('Lag')
            plt.ylabel('Wartość')
            plt.legend()
            plt.grid(True)
            autocorr_png = f"{os.path.splitext(file_name)[0]}_pcal{idx}_autocorr.png"
            autocorr_png_path = os.path.join(output_dir, autocorr_png)
            plt.savefig(autocorr_png_path)
            plt.close()
            logging.debug(f"Wykres funkcji autokorelacyjnej zapisany jako {autocorr_png_path}")

            # Generowanie wykresu widma mocy
            plt.figure(figsize=(10, 6))
            plt.plot(freqs, power_spectrum, label='Widmo Mocy')
            plt.title(f'Widmo Mocy {idx} - {file_name}')
            plt.xlabel('Częstotliwość (Hz)')
            plt.ylabel('Moc (dB)')
            plt.legend()
            plt.grid(True)
            fft_png = f"{os.path.splitext(file_name)[0]}_pcal{idx}_fft.png"
            fft_png_path = os.path.join(output_dir, fft_png)
            plt.savefig(fft_png_path)
            plt.close()
            logging.debug(f"Wykres widma mocy zapisany jako {fft_png_path}")

            # Generowanie interaktywnego wykresu widma mocy za pomocą Plotly
            trace = go.Scatter(
                x=freqs,
                y=power_spectrum,
                mode='lines',
                name='Widmo Mocy'
            )
            layout = go.Layout(
                title=f'Interaktywne Widmo Mocy {idx} - {file_name}',
                xaxis=dict(title='Częstotliwość (Hz)'),
                yaxis=dict(title='Moc (dB)')
            )
            fig = go.Figure(data=[trace], layout=layout)
            interactive_html = f"{os.path.splitext(file_name)[0]}_pcal{idx}_interactive.html"
            interactive_html_path = os.path.join(output_dir, interactive_html)
            fig.write_html(interactive_html_path)
            logging.debug(f"Interaktywny wykres widma mocy zapisany jako {interactive_html_path}")

            # Dodanie informacji o wygenerowanych plikach do listy
            plots_info.append({
                'autocorr_png': autocorr_png,
                'fft_png': fft_png,
                'interactive_html': interactive_html
            })
        except Exception as e:
            logging.error(f"Błąd podczas przetwarzania funkcji autokorelacyjnej {idx} w pliku {file_name}: {e}")
            continue

    return {
        'file_name': file_name,
        'plots': plots_info
    }

def generate_report(results, output_dir, report_file):
    """
    Tworzy raport HTML zawierający wszystkie wyniki analizy.

    Raport zawiera:
    - Nazwę pliku danych.
    - Wykres funkcji autokorelacyjnej (PNG).
    - Wykres widma mocy (PNG).
    - Link do interaktywnego wykresu widma mocy (HTML).

    Args:
        results (list): Lista wyników analizy dla każdego pliku.
        output_dir (str): Katalog, w którym zapisać raport.
        report_file (str): Nazwa pliku raportu.
    """
    report_path = os.path.join(output_dir, report_file)
    try:
        with open(report_path, 'w', encoding='utf-8') as report:
            report.write("<html><head><meta charset='UTF-8'><title>Raport Spectrum Analysis EP</title></head><body>\n")
            report.write("<h1>Raport Spectrum Analysis EP</h1>\n")
            report.write(f"<p>Data generacji raportu: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n")
            for result in results:
                if result is None:
                    continue
                report.write("<div class='file-section'>\n")
                report.write(f"    <h2>Plik danych: {result['file_name']}</h2>\n")
                for i, plot in enumerate(result['plots'], start=1):
                    report.write(f"    <h3>Funkcja autokorelacyjna {i}</h3>\n")
                    report.write(f"    <p><strong>Funkcja autokorelacyjna {i}:</strong></p>\n")
                    report.write(f"    <img src='{plot['autocorr_png']}' alt='Autokorelacja {i}' width='800'>\n")
                    report.write(f"    <p><strong>Widmo mocy {i}:</strong></p>\n")
                    report.write(f"    <img src='{plot['fft_png']}' alt='Widmo mocy {i}' width='800'>\n")
                    report.write(f"    <p><strong>Interaktywny wykres widma mocy {i}:</strong> <a href='{plot['interactive_html']}' target='_blank'>Otwórz</a></p>\n")
                    report.write("    <hr>\n")
                report.write("</div>\n")
            report.write("</body></html>")
        logging.info(f"Raport HTML został zapisany jako {report_path}")
    except Exception as e:
        logging.error(f"Błąd podczas generowania raportu HTML: {e}")

def main():
    """
    Główna funkcja skryptu.
    Wykonuje całą analizę spektralną na podstawie konfiguracji.
    """
    # Parsowanie argumentów wiersza poleceń
    args = parse_arguments()
    config = load_config(args.config)

    # Ustawienie globalnych zmiennych z konfiguracji
    global fs
    fs = config.get('fs', 1000)
    data_dir = config.get('data_dir', 'data')
    exclude_file = config.get('exclude_file', os.path.join(data_dir, 'exclude'))
    output_dir = config.get('output_dir', 'output')
    num_functions = config.get('num_functions', 4)
    samples_per_function = config.get('samples_per_function', 4097)
    log_file = config.get('log_file', 'spectrum_analysis.log')
    log_level = config.get('log_level', 'INFO')
    fft_type = config.get('fft_type', 'rfft')
    report_file = config.get('report_file', 'raport.html')

    # Konfiguracja logowania
    setup_logging(log_file, log_level)

    # Logowanie informacji startowych
    logging.info("Rozpoczynanie analizy Spectrum Analysis EP")
    logging.info(f"Katalog z danymi wejściowymi: {data_dir}")
    logging.info(f"Plik wykluczeń: {exclude_file}")
    logging.info(f"Katalog z wynikami: {output_dir}")
    logging.info(f"Częstotliwość próbkowania: {fs} Hz")
    logging.info(f"Liczba funkcji autokorelacyjnych: {num_functions}")
    logging.info(f"Liczba próbek na funkcję: {samples_per_function}")
    logging.info(f"Typ FFT: {fft_type}")

    # Utworzenie folderu output, jeśli nie istnieje
    os.makedirs(output_dir, exist_ok=True)

    # Wczytanie listy plików do wykluczenia
    exclude_files = load_exclude_list(exclude_file)

    # Pobranie listy plików do analizy
    data_files = get_data_files(data_dir, exclude_files)

    if not data_files:
        logging.warning("Brak plików do analizy. Skrypt zostanie zakończony.")
        return

    # Przygotowanie argumentów do przetwarzania równoległego
    process_args = [
        (file_path, output_dir, fs, num_functions, samples_per_function, fft_type)
        for file_path in data_files
    ]

    # Przetwarzanie plików równolegle przy użyciu puli procesów
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(generate_plots, process_args)

    # Filtracja wyników (usunięcie None)
    results = [result for result in results if result is not None]

    if not results:
        logging.warning("Brak wyników do raportu. Skrypt zostanie zakończony.")
        return

    # Generowanie raportu HTML
    generate_report(results, output_dir, report_file)
    logging.info("Analiza zakończona.")

if __name__ == "__main__":
    main()
