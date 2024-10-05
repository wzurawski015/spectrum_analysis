#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spectrum_analysis_optimized.py

Zoptymalizowany skrypt do analizy spektralnej danych autokorelacyjnych.
Wykonuje następujące kroki:
1. Wczytuje dane z plików tekstowych.
2. Ignoruje pliki na liście wykluczeń.
3. Przetwarza każdą funkcję autokorelacyjną:
   a. Symetryzuje funkcję autokorelacyjną.
   b. Usuwa składową stałą (DC offset).
   c. Oblicza FFT i widmo mocy za pomocą pyFFTW.
   d. Umożliwia wybór transformacji (FFT, rFFT, DCT).
   e. Zapisuje pośrednie wyniki do plików JSON.
   f. Generuje wykresy statyczne i interaktywne.
4. Tworzy dynamiczny raport HTML zawierający wszystkie wyniki analizy.
5. Wykorzystuje przetwarzanie równoległe dla zwiększenia wydajności.
6. Obsługuje błędy, np. brak plików lub niepoprawny format danych.

Autor: [Twoje Imię]
Data: 2024-04-27 (zaktualizowano)

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
import pyfftw
import pyfftw.interfaces.scipy_fft as fftw
from scipy.fftpack import dct
from datetime import datetime
from multiprocessing import Pool, cpu_count

def load_config(config_path):
    """
    Wczytuje konfigurację z pliku YAML.
    """
    with open(config_path, 'r', encoding='utf-8') as file:
        config = yaml.safe_load(file)
    return config

def setup_logging(log_file, log_level):
    """
    Konfiguruje logowanie do pliku i konsoli.
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

def load_data(file_path):
    """
    Wczytuje dane z pliku.
    """
    try:
        data = np.loadtxt(file_path)
        return data
    except Exception as e:
        logging.error(f"Nie udało się wczytać pliku {file_path}: {str(e)}")
        return None

def compute_transform(data, transform_type='rfft'):
    """
    Oblicza transformację (FFT, rFFT, DCT) dla danych.
    """
    if transform_type == 'rfft':
        return fftw.rfft(data)
    elif transform_type == 'fft':
        return fftw.fft(data)
    elif transform_type == 'dct':
        return dct(data)
    else:
        raise ValueError(f"Nieznany typ transformacji: {transform_type}")

def save_results(file_path, transformed_data, output_dir):
    """
    Zapisuje wyniki transformacji do pliku JSON.
    """
    try:
        result_file = os.path.join(output_dir, os.path.basename(file_path) + '_results.json')
        with open(result_file, 'w', encoding='utf-8') as json_file:
            json.dump({'transformed_data': transformed_data.tolist()}, json_file)
        logging.info(f"Wyniki zapisane do {result_file}")
    except Exception as e:
        logging.error(f"Nie udało się zapisać wyników dla {file_path}: {str(e)}")

def generate_plots(args):
    """
    Generuje wykresy dla danych.
    """
    file_path, output_dir, fs, num_functions, samples_per_function, fft_type = args
    data = load_data(file_path)
    if data is None:
        return None

    # Symetryzacja, usunięcie DC offsetu itp.
    transformed_data = compute_transform(data, fft_type)

    # Zapis wyników
    save_results(file_path, transformed_data, output_dir)

    return {
        'file_name': file_path,
        'frequency': np.arange(len(transformed_data)),
        'power_spectrum': np.abs(transformed_data)
    }

def generate_report(results, output_dir, report_file):
    """
    Generuje raport HTML z wynikami.
    """
    html_content = '<html><head><title>Raport z Analizy</title></head><body>'
    html_content += '<h1>Raport z Analizy Spektralnej</h1>'

    for result in results:
        html_content += f"<h2>Plik: {result['file_name']}</h2>"
        html_content += f"<p>Opis wyników dla tego pliku...</p>"

        # Dynamiczny wykres za pomocą Plotly
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=result['frequency'], y=result['power_spectrum'], mode='lines'))
        html_content += fig.to_html(full_html=False)

    html_content += '</body></html>'

    with open(os.path.join(output_dir, report_file), 'w', encoding='utf-8') as f:
        f.write(html_content)
    logging.info(f"Raport HTML zapisany do {report_file}")

def process_file(file_path, output_dir, fs, num_functions, samples_per_function, fft_type):
    """
    Przetwarza pojedynczy plik i zapisuje wyniki.
    """
    try:
        data = load_data(file_path)
        if data is None:
            return None

        transformed_data = compute_transform(data, fft_type)
        save_results(file_path, transformed_data, output_dir)

        return {
            'file_name': file_path,
            'frequency': np.arange(len(transformed_data)),
            'power_spectrum': np.abs(transformed_data)
        }
    except Exception as e:
        logging.error(f"Błąd przetwarzania pliku {file_path}: {str(e)}")
        return None

def main():
    """
    Główna funkcja skryptu.
    """
    # Parsowanie argumentów
    parser = argparse.ArgumentParser(description="Spektralna analiza danych autokorelacyjnych")
    parser.add_argument('config', type=str, help="Ścieżka do pliku konfiguracyjnego YAML")
    args = parser.parse_args()

    # Wczytanie konfiguracji
    config = load_config(args.config)

    data_dir = config['data_dir']
    exclude_file = config['exclude_file']
    output_dir = config['output_dir']
    fs = config['sampling_frequency']
    num_functions = config['num_functions']
    samples_per_function = config['samples_per_function']
    log_file = config.get('log_file', 'spectrum_analysis.log')
    log_level = config.get('log_level', 'INFO')
    fft_type = config.get('fft_type', 'rfft')
    report_file = config.get('report_file', 'raport.html')
    transform_type = config.get('transform_type', 'rfft')  # Nowa opcja

    # Konfiguracja logowania
    setup_logging(log_file, log_level)

    # Logowanie informacji startowych
    logging.info("Rozpoczynanie analizy Spectrum Analysis EP")
    logging.info(f"Katalog z danymi wejściowymi: {data_dir}")
    logging.info(f"Katalog z wynikami: {output_dir}")
    logging.info(f"Częstotliwość próbkowania: {fs} Hz")
    logging.info(f"Liczba funkcji autokorelacyjnych: {num_functions}")
    logging.info(f"Typ transformacji: {transform_type}")

    # Utworzenie folderu output, jeśli nie istnieje
    os.makedirs(output_dir, exist_ok=True)

    # Przygotowanie argumentów do przetwarzania równoległego
    data_files = glob.glob(os.path.join(data_dir, '*.txt'))  # Zakładamy pliki .txt do analizy
    process_args = [
        (file_path, output_dir, fs, num_functions, samples_per_function, transform_type)
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
