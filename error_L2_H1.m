function [errL2, errH1] = error_L2_H1(u, coordinates, elements3)

% ERROR_L2_H1  Compute L2 and H1 errors via 3‑point Gaussian on each triangle
%
%   [errL2,errH1] = error_L2_H1(u, coords, elems3)
%   u          : FEM nodal values (Nnodes×1)
%   coordinates: Nnodes×2 array of (x,y)
%   elements3  : Nelems×3 connectivity (node indices)

 
                

w7 = [1/6; 1/6; 1/6];           
b7 =   [2/3, 1/6;   % (ξ,η)
         1/6, 2/3;
         1/6, 1/6];

 L2sq  = 0;
  H1sq  = 0;

  for k = 1:size(elements3,1)
  
  % vertex coordinates of K
  verts = coordinates(elements3(k,:), :);   % 3×2 array
  uK    = u(elements3(k,:));                % 3×1 nodal values

  % affine map: x = a1 + B*[ξ;η]
  a1   = verts(1,:)';                       % column vector [x1; y1]
  J    = [ (verts(2,:)-verts(1,:))', ...     % 2×2 matrix
           (verts(3,:)-verts(1,:))' ];
  detJ = abs(det(J));

  % ∇u_h on K is constant: ∇u_h = B^{-T} * [u2−u1; u3−u1]
  grad_uh = J' \ ([uK(2)-uK(1); uK(3)-uK(1)])


    % loop over quad points
    for q = 1:3
      lam1 = b7(q,1);
      lam2 = b7(q,2);
      %
        %a1   = verts(1,:)';   
      % Physical quadrature point
      xq = a1 + J*[lam1; lam2];    
      % exact solution and its gradient
      u_ex   =  sin(pi*xq(1)).*sin(pi*xq(2));
      grad_ex =[pi*cos(pi*xq(1))*sin(pi*xq(2)), ...
                      pi*sin(pi*xq(1))*cos(pi*xq(2))];

      % FEM solution at xq: uh_q = [lam3 lam1 lam2] * u(v)
      uh_q =  uK(1)*(1-lam1-lam2) + uK(2)*lam1 + uK(3)*lam2;

      

        % accumulate L2 error
        L2sq = L2sq + detJ * w7(q) * (u_ex - uh_q)^2;
    
    
        % accumulate H1-semi error
        diffg = grad_ex' - grad_uh;
        H1sq  = H1sq  + detJ * w7(q) * (diffg' * diffg);
    end

  errL2 = sqrt(L2sq);
  errH1 = sqrt(H1sq);
end
