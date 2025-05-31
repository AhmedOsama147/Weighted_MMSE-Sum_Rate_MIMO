% Recursive function for sphere decoding
function x_hat = search_level(n, dist, d, ref_ang, x_partial, best_distance, R)
for i = n:-1:1
    for m = 1:length(ref_ang)
        x_partial(i) = ref_ang(m);
        dist = (d(i) - R(i,i:n)*x_partial(i:n))^2 + dist;
        if dist <= best_distance
            best_distance = dist;
            break;
        end
    end
end
x_hat = x_partial;
end