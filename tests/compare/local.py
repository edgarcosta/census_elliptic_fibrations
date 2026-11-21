def fibres(text):
    if text == "[]":
        return []
    rows = []
    for entry in text.split(";"):
        place, degree, kind, index, vd, conductor, c, trace, cycles = entry.split(":")
        rows.append((place, int(degree), kind, int(index), int(vd),
                     int(conductor), int(c), int(trace),
                     [] if cycles == "[]" else list(map(int, cycles.split()))))
    return rows


def check_identities(text, height):
    rows = fibres(text)
    assert sum(r[1] * r[4] for r in rows) == 12 * height, rows
    for place, d, kind, n, vd, conductor, c, trace, cycles in rows:
        components = n if kind == "I" else n + 5 if kind == "I*" else {
            "II": 1, "III": 2, "IV": 3, "IV*": 7, "III*": 8, "II*": 9,
        }[kind]
        assert conductor == vd + 1 - components, (place, kind, conductor, vd)
        assert sum(cycles) == d * (components - 1), (place, kind, cycles)
        assert cycles == sorted(cycles) and all(k % d == 0 for k in cycles)
