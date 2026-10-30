import json, sys
m = json.load(open(sys.argv[1]))
for comp, d in m["noninferiority"].items():
    for lv in ("peptidoform", "protein"):
        x = d[lv]
        s = " ".join(f"own{a}={(x['own_q'][a]['ratio'] or float('nan')):.3f}[{x['own_q'][a]['ci95'][0]:.3f},{x['own_q'][a]['ci95'][1]:.3f}]({x['own_q'][a]['comparator_verdict_deciding']})" for a in ("0.001","0.01","0.05"))
        s2 = " ".join(f"m{f}={(x['matched_fdp'][f]['ratio'] or float('nan')):.3f}[{x['matched_fdp'][f]['ci95'][0]:.3f},{x['matched_fdp'][f]['ci95'][1]:.3f}]" for f in sorted(x["matched_fdp"]))
        print(comp, lv, s2, "|", s)
for meth, r in m["methods"].items():
    if "verdicts" not in r: continue
    v = r["verdicts"]["peptidoform"]
    print(meth, {a: (round(v[a]["P"]["a"]["mean_fdp"], 4), round(v[a]["C"]["a"]["mean_fdp"], 4), v[a]["deciding"]) for a in ("0.001", "0.01", "0.05")},
          {a: round(r["ids_dataset_level"]["peptidoform"][a]["mean"]) for a in ("0.001","0.01","0.05")})
