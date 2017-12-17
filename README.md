# tsg2_zal_2
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

```bash
virtualenv -p /usr/bin/python3.6 ~/venvs/tsg_2_2
source ~/venvs/tsg_2_2/bin/activate
pip install --upgrade pip
pip install pip-tools
pip-sync
```

## Tesy
TODO

## Uruchomienie programu
`./assembly input.fasta output.fasta`
