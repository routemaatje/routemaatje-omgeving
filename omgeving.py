#!/usr/bin/env python3
"""Routemaatje: omgeving van campings berekenen uit OpenStreetMap-gegevens (Geofabrik).

Gebruik:
    omgeving.py LAND export.geojsonseq [routes_wandel.geojsonseq routes_fiets.geojsonseq] -o uit.ndjson.gz [--datum JJJJMMDD]

- export.geojsonseq: uitvoer van `osmium export` (punten, lijnen en vlakken met tags en id "type+nummer").
- routes_*: lijnen van de leden van wandel- en fietsroutes (alle lijnen in die bestanden horen bij een route).

Uitvoer (gzip, een regel per camping of camperplaats):
  eerste regel : {"_meta": {"land", "aantal", "bron_datum", "versie", "gemaakt"}}
  daarna       : {"t": "n|w|r", "i": osm-id, "lat": .., "lon": .., "soort": "camping|camperplaats",
                  "tellingen": {"water": {"250": n, "1000": n, "3000": n}, "rivier": {...}, ...}}

Een telling is het aantal verschillende objecten waarvan een deel binnen die afstand (in meters) ligt.
Ligt de camping IN een vlak (meer, bos, natuurgebied), dan telt dat object op elke afstand mee.
Alleen standaardbibliotheek, numpy en scipy.
"""
import argparse
import gzip
import json
import math
import sys
from datetime import datetime, timezone

import numpy as np
from scipy.spatial import cKDTree

VERSIE = 2
AARDE = 6371008.8

# categorie -> afstanden in meter
AFSTANDEN = {
    "water": [250, 1000, 3000],
    "rivier": [250, 1000, 3000],
    "zee": [250, 1000, 3000],
    "bos": [250, 1000, 3000],
    "natuur": [250, 1000, 3000],
    "berg": [10000],
    "supermarkt": [1000, 3000],
    "restaurant": [1000],
    "wandel": [2000],
    "fiets": [2000],
}

GEEN_MEER = {"river", "stream", "canal", "ditch", "drain", "basin", "pond", "wastewater",
             "swimming_pool", "fountain", "moat", "stream_pool", "lock", "reflecting_pool"}


def categorieen(t, geom):
    """Welke categorieen past dit object? t = tags, geom = 'Point', 'LineString' of 'Polygon'."""
    uit = []
    nat = t.get("natural")
    if geom in ("Polygon", "MultiPolygon", "Point") and (
        (nat == "water" and t.get("water") not in GEEN_MEER) or t.get("landuse") == "reservoir"
    ):
        uit.append("water")
    if t.get("waterway") == "river" or (nat == "water" and t.get("water") == "river"):
        uit.append("rivier")
    if nat == "coastline" or (nat == "beach" and t.get("surface") != "artificial"):
        uit.append("zee")
    if geom in ("Polygon", "MultiPolygon") and (nat == "wood" or t.get("landuse") == "forest"):
        uit.append("bos")
    if geom in ("Polygon", "MultiPolygon", "Point") and (
        t.get("leisure") == "nature_reserve"
        or t.get("boundary") == "national_park"
        or (t.get("boundary") == "protected_area" and str(t.get("protect_class")) in {"1", "2", "3", "4"})
    ):
        uit.append("natuur")
    if nat in ("peak", "volcano") and geom == "Point":
        uit.append("berg")
    if t.get("shop") == "supermarket":
        uit.append("supermarkt")
    if t.get("amenity") == "restaurant":
        uit.append("restaurant")
    return uit


def xyz(lat, lon):
    """Breedte/lengte (graden) naar punten op een bol; afstanden zijn dan gewone rechte afstanden (m)."""
    la = np.radians(np.asarray(lat, dtype=float))
    lo = np.radians(np.asarray(lon, dtype=float))
    c = np.cos(la)
    return np.column_stack((AARDE * c * np.cos(lo), AARDE * c * np.sin(lo), AARDE * np.sin(la)))


def hav(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    a = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(math.radians(lon2 - lon1) / 2) ** 2
    return 2 * AARDE * math.asin(math.sqrt(a))


def punten_langs(coords, max_stap=250.0):
    """Hoekpunten van een lijn of ring, aangevuld zodat de punten hoogstens ongeveer max_stap meter uit elkaar liggen."""
    uit = []
    vorige = None
    for c in coords:
        lon, lat = c[0], c[1]
        if vorige is not None:
            d = hav(vorige[1], vorige[0], lat, lon)
            if d > max_stap:
                n = int(d // max_stap)
                for k in range(1, n + 1):
                    f = k / (n + 1)
                    uit.append((vorige[1] + (lat - vorige[1]) * f, vorige[0] + (lon - vorige[0]) * f))
        uit.append((lat, lon))
        vorige = (lon, lat)
    return uit


def in_ring(lon, lat, ring):
    """Punt in ring (ray casting)."""
    binnen = False
    j = len(ring) - 1
    for i in range(len(ring)):
        xi, yi = ring[i][0], ring[i][1]
        xj, yj = ring[j][0], ring[j][1]
        if (yi > lat) != (yj > lat) and lon < (xj - xi) * (lat - yi) / (yj - yi) + xi:
            binnen = not binnen
        j = i
    return binnen


def in_vlak(lon, lat, ringen):
    """ringen[0] = buitenrand, rest = gaten."""
    if not ringen or not in_ring(lon, lat, ringen[0]):
        return False
    return not any(in_ring(lon, lat, g) for g in ringen[1:])


def lees_geojsonseq(pad):
    opener = gzip.open if str(pad).endswith(".gz") else open
    with opener(pad, "rt", encoding="utf-8") as f:
        for regel in f:
            regel = regel.strip().lstrip("\x1e").strip()
            if not regel:
                continue
            try:
                yield json.loads(regel)
            except ValueError:
                continue


def polygonen(geom):
    if geom["type"] == "Polygon":
        return [geom["coordinates"]]
    if geom["type"] == "MultiPolygon":
        return geom["coordinates"]
    return []


def midden(ringen):
    ring = ringen[0]
    lons = [c[0] for c in ring]
    lats = [c[1] for c in ring]
    return (min(lats) + max(lats)) / 2, (min(lons) + max(lons)) / 2


def split_id(fid):
    """osmium-id 'w123' / 'n5' / 'r9' -> ('w', 123)."""
    s = str(fid)
    if s and s[0] in "nwra" and s[1:].isdigit():
        n = int(s[1:])
        if s[0] == "a":   # osmium-vlak: even = way (nummer / 2), oneven = relatie ((nummer - 1) / 2)
            return ("w", n // 2) if n % 2 == 0 else ("r", n // 2)
        return s[0], n
    return None, None


class Verzamelaar:
    """Houdt per categorie de steekpunten (met object-nummer) en de vlakken bij."""

    def __init__(self):
        self.pt = {c: ([], [], []) for c in AFSTANDEN}   # lat, lon, object-nr
        self.vlakken = {c: [] for c in ("water", "bos", "natuur")}   # (nr, bbox, ringen)
        self.nr = 0
        self.campings = []   # (t, i, lat, lon, soort)

    def voeg_punten(self, cat, punten):
        self.nr += 1
        la, lo, nr = self.pt[cat]
        for lat, lon in punten:
            la.append(lat)
            lo.append(lon)
            nr.append(self.nr)

    def verwerk(self, ft):
        geom = ft.get("geometry")
        if not geom:
            return
        t = ft.get("properties") or {}
        fid = ft.get("id", (ft.get("properties") or {}).get("@id"))
        gt = geom["type"]
        soort = None
        if t.get("tourism") == "camp_site":
            soort = "camping"
        elif t.get("tourism") == "caravan_site":
            soort = "camperplaats"
        if soort and gt in ("Point", "Polygon", "MultiPolygon"):
            if gt == "Point":
                lon, lat = geom["coordinates"][0], geom["coordinates"][1]
            else:
                pol = polygonen(geom)
                lat, lon = midden(pol[0])
            ty, nm = split_id(fid)
            if ty:
                self.campings.append((ty, nm, lat, lon, soort))
        for cat in categorieen(t, "Polygon" if gt == "MultiPolygon" else gt):
            if gt == "Point":
                self.voeg_punten(cat, [(geom["coordinates"][1], geom["coordinates"][0])])
            elif gt == "LineString":
                self.voeg_punten(cat, punten_langs(geom["coordinates"]))
            elif gt in ("Polygon", "MultiPolygon"):
                for ringen in polygonen(geom):
                    self.voeg_punten(cat, punten_langs(ringen[0]))
                    if cat in self.vlakken:
                        lons = [c[0] for c in ringen[0]]
                        lats = [c[1] for c in ringen[0]]
                        self.vlakken[cat].append((self.nr, (min(lats), min(lons), max(lats), max(lons)), ringen))

    def voeg_route_lijnen(self, cat, pad):
        for ft in lees_geojsonseq(pad):
            geom = ft.get("geometry") or {}
            if geom.get("type") == "LineString":
                self.voeg_punten(cat, punten_langs(geom["coordinates"]))


def bereken(v):
    """Geeft per camping de tellingen."""
    n = len(v.campings)
    tel = [dict() for _ in range(n)]
    if n == 0:
        return tel
    cxyz = xyz([c[2] for c in v.campings], [c[3] for c in v.campings])
    for cat, radii in AFSTANDEN.items():
        la, lo, nr = v.pt[cat]
        sets = [[set() for _ in radii] for _ in range(n)]
        if la:
            boom = cKDTree(xyz(la, lo))
            nr_arr = np.asarray(nr)
            ruim = boom.query_ball_point(cxyz, r=max(radii))
            dist_cache = None
            for k in range(n):
                idx = ruim[k]
                if not idx:
                    continue
                idx = np.asarray(idx)
                d = np.linalg.norm(boom.data[idx] - cxyz[k], axis=1)
                for ri, r in enumerate(radii):
                    sets[k][ri].update(nr_arr[idx[d <= r]].tolist())
        # camping binnen een vlak: telt op alle afstanden
        if cat in v.vlakken and v.vlakken[cat]:
            clat = np.asarray([c[2] for c in v.campings])
            clon = np.asarray([c[3] for c in v.campings])
            kb = cKDTree(np.column_stack((clat, clon)))
            for (vnr, (mi_la, mi_lo, ma_la, ma_lo), ringen) in v.vlakken[cat]:
                mid = ((mi_la + ma_la) / 2, (mi_lo + ma_lo) / 2)
                straal = math.hypot(ma_la - mi_la, ma_lo - mi_lo) / 2 + 1e-9
                for k in kb.query_ball_point(mid, r=straal):
                    if mi_la <= clat[k] <= ma_la and mi_lo <= clon[k] <= ma_lo and in_vlak(clon[k], clat[k], ringen):
                        for ri in range(len(radii)):
                            sets[k][ri].add(vnr)
        for k in range(n):
            tel[k][cat] = {str(r): len(sets[k][ri]) for ri, r in enumerate(radii)}
    return tel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("land")
    ap.add_argument("export")
    ap.add_argument("routes_wandel", nargs="?")
    ap.add_argument("routes_fiets", nargs="?")
    ap.add_argument("-o", "--uit", required=True)
    ap.add_argument("--datum", default=datetime.now(timezone.utc).strftime("%Y%m%d"))
    ap.add_argument("--kern", default=None, help="minlon,minlat,maxlon,maxlat: alleen campings met middelpunt in deze doos (tegels)")
    a = ap.parse_args()

    v = Verzamelaar()
    for ft in lees_geojsonseq(a.export):
        v.verwerk(ft)
    if a.routes_wandel:
        v.voeg_route_lijnen("wandel", a.routes_wandel)
    if a.routes_fiets:
        v.voeg_route_lijnen("fiets", a.routes_fiets)

    if a.kern:
        x0, y0, x1, y1 = [float(s) for s in a.kern.split(",")]
        v.campings = [c for c in v.campings if x0 <= c[3] < x1 and y0 <= c[2] < y1]

    tel = bereken(v)
    with gzip.open(a.uit, "wt", encoding="utf-8") as f:
        meta = {"land": a.land.upper(), "aantal": len(v.campings), "bron_datum": a.datum, "versie": VERSIE,
                "gemaakt": datetime.now(timezone.utc).isoformat(timespec="seconds"), "afstanden": AFSTANDEN}
        f.write(json.dumps({"_meta": meta}, ensure_ascii=False) + "\n")
        for (ty, nm, lat, lon, soort), t in zip(v.campings, tel):
            f.write(json.dumps({"t": ty, "i": nm, "lat": round(lat, 6), "lon": round(lon, 6), "soort": soort,
                                "tellingen": t}, ensure_ascii=False, separators=(",", ":")) + "\n")
    print(f"{a.land.upper()}: {len(v.campings)} campings/camperplaatsen berekend", file=sys.stderr)


if __name__ == "__main__":
    main()
