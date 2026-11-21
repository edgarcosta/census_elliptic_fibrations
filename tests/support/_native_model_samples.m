function model_samples(R)
    s := R.1;
    t := R.2;
    p := Characteristic(R);
    if p eq 2 then
        samples := [
            [R | s*t, 0, s^6+t^6, 0, 0],
            [R | s^2, 0, 0, 0, s^12],
            [R | 0, 0, 0, s^7*t+s*t^7, 0],
            [R | 0, 0, 0, 0, s*t^11],
            [R | s*t, 0, s^5*t+s*t^5, s^8+t^8, 0],
            [R | 0, 0, 0, 0, 0]];
    else
        samples := [
            [R | 0, s^3*t+s*t^3, 0, s^8+t^8, 0],
            [R | 0, s^4, 0, 0, s^12],
            [R | 0, 0, 0, 0, s^11*t+s*t^11],
            [R | 0, 0, 0, 0, s*t^11],
            [R | 0, s^3*t+s*t^3, 0, 0, 0],
            [R | 0, 0, 0, 0, 0]];
    end if;
    return samples;
end function;
