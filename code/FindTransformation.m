%function of calculation of A and B
% I/P
% Fp = Column of 5 facial features (64*64)
% F_cordinate = Row of 5 features cordinate 
% o/p = get A and B

function [A,B] = FindTransformation(Fp, F_cordinate)
    def_matrix = [F_cordinate(1:2) 0 0 1 0; 0 0 F_cordinate(1:2) 0 1;
        F_cordinate(3:4) 0 0 1 0; 0 0 F_cordinate(3:4) 0 1;
        F_cordinate(5:6) 0 0 1 0; 0 0 F_cordinate(5:6) 0 1;
        F_cordinate(7:8) 0 0 1 0; 0 0 F_cordinate(7:8) 0 1;
        F_cordinate(9:10) 0 0 1 0; 0 0 F_cordinate(9:10) 0 1;];
    x = def_matrix\Fp; %Least square solution of the system created
    A = [x(1) x(2); x(3) x(4)]; % Defining the A and b transformation
    B = [x(5); x(6)];
end