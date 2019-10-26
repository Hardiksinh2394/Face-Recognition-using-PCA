
function F_prime = ApplyTransformation(A,B,F_cordinate) % apply a and b to get f prime 
    F_x = F_cordinate(1:2:end);
    F_y = F_cordinate(2:2:end);
    F_prime = [];
    for i = 1:5
        f_feature = [F_x(i) F_y(i)]';
        fi_feature = A*f_feature+B;
        F_prime = [F_prime fi_feature'];
    end
end