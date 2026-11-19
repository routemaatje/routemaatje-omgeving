import gzip, json, subprocess, sys, os, tempfile, math
d = tempfile.mkdtemp()
def deg_n(m): return m / 111195.0           # meters noord -> graden breedte
def deg_e(m, lat=52.0): return m / (111195.0 * math.cos(math.radians(lat)))
lat0, lon0 = 52.0, 5.0
def F(id_, tags, geom): return {"type": "Feature", "id": id_, "properties": tags, "geometry": geom}
def pt(n, e): return {"type": "Point", "coordinates": [lon0 + deg_e(e), lat0 + deg_n(n)]}
def box(n1, e1, n2, e2): return {"type": "Polygon", "coordinates": [[[lon0+deg_e(e1), lat0+deg_n(n1)], [lon0+deg_e(e2), lat0+deg_n(n1)], [lon0+deg_e(e2), lat0+deg_n(n2)], [lon0+deg_e(e1), lat0+deg_n(n2)], [lon0+deg_e(e1), lat0+deg_n(n1)]]]}
def line(*p): return {"type": "LineString", "coordinates": [[lon0+deg_e(e), lat0+deg_n(n)] for n, e in p]}
feats = [
 F("n1", {"tourism": "camp_site", "name": "Test"}, pt(0, 0)),
 F("w2", {"natural": "water", "water": "lake"}, box(-100, -100, 100, 100)),        # camping ligt in het meer
 F("w3", {"waterway": "river"}, line((-2000, 500), (2000, 500))),                      # rivier op 500 m
 F("w4", {"landuse": "forest"}, box(1900, -300, 2400, 300)),                           # bos op ~1,9 km
 F("w5", {"natural": "coastline"}, line((-5000, 2800), (5000, 2800))),                 # kust op 2,8 km
 F("n6", {"shop": "supermarket"}, pt(800, 0)),
 F("n7", {"amenity": "restaurant"}, pt(300, 300)),
 F("n8", {"amenity": "restaurant"}, pt(5000, 0)),                                      # te ver
 F("n9", {"natural": "peak"}, pt(0, 7000)),
 F("w10", {"leisure": "nature_reserve"}, box(-5000, -5000, 5000, -150)),               # natuurgebied grenst op 150 m
 F("w11", {"natural": "water", "water": "pond"}, box(-50, 900, 50, 950)),              # vijver telt niet als meer
]
p = os.path.join(d, "export.geojsonseq")
open(p, "w").write("".join("\x1e" + json.dumps(f) + "\n" for f in feats))
r = os.path.join(d, "route.geojsonseq")
open(r, "w").write("\x1e" + json.dumps(F("w20", {}, line((-1000, 1500), (1000, 1500)))) + "\n"
                   + "\x1e" + json.dumps(F("w21", {}, line((-1000, 3000), (1000, 3000)))) + "\n")
uit = os.path.join(d, "uit.ndjson.gz")
subprocess.check_call([sys.executable, os.path.join(os.path.dirname(__file__), "..", "omgeving.py"), "NL", p, r, r, "-o", uit, "--datum", "20261119"])
regels = [json.loads(x) for x in gzip.open(uit, "rt")]
meta, c = regels[0]["_meta"], regels[1]
assert meta["aantal"] == 1 and meta["land"] == "NL", meta
t = c["tellingen"]
assert (c["t"], c["i"], c["soort"]) == ("n", 1, "camping"), c
assert t["water"] == {"250": 1, "1000": 1, "3000": 1}, t["water"]          # alleen het meer, vijver niet
assert t["rivier"] == {"250": 0, "1000": 1, "3000": 1}, t["rivier"]
assert t["bos"] == {"250": 0, "1000": 0, "3000": 1}, t["bos"]
assert t["zee"] == {"250": 0, "1000": 0, "3000": 1}, t["zee"]
assert t["supermarkt"] == {"1000": 1, "3000": 1}, t["supermarkt"]
assert t["restaurant"] == {"1000": 1}, t["restaurant"]
assert t["berg"] == {"10000": 1}, t["berg"]
assert t["natuur"] == {"250": 1, "1000": 1, "3000": 1}, t["natuur"]
assert t["wandel"] == {"2000": 1}, t["wandel"]                               # alleen de lijn op 1,5 km
print("ALLE TESTS OK", json.dumps(t))
