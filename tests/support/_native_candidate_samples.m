function candidate_samples(F)
    t := F.1;
    if Characteristic(F) eq 2 then
        return [
            <[F | 0, 0, 1, t^6, 0], [1]>,
            <[F | 0, 0, 1, t^7, 0], [1, 4]>,
            <[F | t, 0, t^6+1, 0, 0], [1]>,
            <[F | t^2, 0, 1, 0, 0], [1]>];
    end if;
    return [
        <[F | 0, 0, 0, t^8, t+1], [1]>,
        <[F | 0, t^3+t, 0, t^8+1, 0], [1, 4]>,
        <[F | 0, 0, 0, 1, t^7], [1, 4, 9]>,
        <[F | 0, t^4+1, 0, t^4, 0], [1]>];
end function;
