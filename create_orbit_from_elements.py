def create_orbit_from_elements(row):
    a = row["半长轴（km)"] * u.km
    ecc = row["偏心率"] * u.one
    inc = row["轨道倾角(°)"] * u.deg
    raan = row["升交点赤经(°)"] * u.deg
    argp = row["近地点幅角(°）"] * u.deg
    nu = row["平近点角(°)"] * u.deg

    return Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)