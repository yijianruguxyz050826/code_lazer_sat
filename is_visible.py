def is_visible(radar_idx, target_idx, time_slot):
    orbit = orbits[target_idx]

    # 使用时间偏移量（相对时间）
    time_of_flight = time_slot * dt_per_minute  # dt_per_minute = 60 * u.s

    try:
        r_target = orbit.propagate(time_of_flight).r.to(u.km)
    except Exception as e:
        print(f"轨道传播失败 at slot {time_slot}: {e}")
        return False

    r_radar = radar_positions[radar_idx]
    distance = np.linalg.norm((r_target - r_radar).value)

    max_distance = sensor_data.loc[radar_idx, "探测距离(km)"]
    return distance <= max_distance