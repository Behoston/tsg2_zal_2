# \[GAS] Genome Assembly as Service

Web service dla projektu Genome Assembly.

## Wymagania
Ten projekt wymaga Pythona w wersji 3.6

Wymagane dodatkowe pakiety:

`pip install -r ./requirements.txt`

## Baza danych
Aby setupować bazę danych SQLite użyj:
`python model.py`

## Uruchomienie aplikacji
Uruchom `python app.py`


## Plany
1. Użyć Dockera
2. Użyć Postgresa zamiast SQLite
3. Używać workerów spełniających założenia co do zasobów (0.5 GB RAM, 1 rdzeń)
