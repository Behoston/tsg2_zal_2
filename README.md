# tsg2_zal_2

[![Build Status](https://travis-ci.org/Behoston/tsg2_zal_2.svg?branch=master)](https://travis-ci.org/Behoston/tsg2_zal_2)

Projekt polega na opracowaniu i zaimplementowaniu algorytmu asemblacji odczytów z sekwencjonowania.

[Opis projektu](https://moodle.mimuw.edu.pl/pluginfile.php?file=%2F24584%2Fmod_resource%2Fcontent%2F0%2Fassignment2.pdf)

## Ograniczenia projektu

- Nie można używać programów ani bibliotek do asemblacji, mapowania, znajdowania uliniowień itp.
- Program musi domyślnie działać na 1 wątku

## Wymagania wydajnościowe
Dla 1000 odczytów o typowych parametrach:
- czas działania nie dłuższy niż 1 h
- zużycie pamięci nie większe niż 0.5 GB

## Środowisko deweloperkskie

- Projekt używa Pythona w wersji 3.6
- Do instalacji pakietów używany jest `pip-sync`
- Wymagany jest `git lfs`
- Wymagany jest `bowtie2` dostępny bezpośrednio z terminala (dodany do PATH)


```bash
virtualenv -p /usr/bin/python3.6 ~/venvs/tsg_2_2
source ~/venvs/tsg_2_2/bin/activate
pip install --upgrade pip
pip install pip-tools
pip-sync
```

## Analiza statyczna
`flake8`

## Zasady pisania kodu
1. Każdy merge request (idealnie każdy commit) musi być zgodny z flake8
2. Każdy merge request powinien zawierać testy jeśli to możliwe/ma sens


## Zasady używania gita
1. Każdy issue na innym branchu
2. Każdy branch jest domyślnie dla jednej osoby (chyba, że wyraźnie jest prośba o współpracę)
3. Merge request musi być zawsze sprawdzony przez przynajmniej 1 osobę (inną niż autor lub autorzy)
4. Merge reuqest musi zawsze odchodzić z mastera (robimy rebase) żeby nie powodować konfliktów

## Testy
Wyniki można porównywać za pomocą skryptu [evaluate_bowtie.sh](./scripts/evaluate_bowtie.sh)

(wymagany jest bowtie2)



```bash

./scripts/evaluate_bowite.sh output.fasta
```

## Uruchomienie programu
`python3.6 assembly.py input.fasta output.fasta`
