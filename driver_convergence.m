clear;  close all;  clc;

ns     = [8 16 32 64,128];
errL2  = zeros(size(ns));
errH1  = zeros(size(ns));

for k = 1:numel(ns)

  % (re)generate the mesh files for this n:
   write_coordinates_dat(ns(k))     
  % plot on log‑log scale
% h = 1./ns(k);
  % solve the Poisson problem on that mesh:
  fem2d  

 
  [errL2, errH1] = error_L2_H1(full(u), coordinates, elements3);

 L2err(k) = errL2;
 H1err(k) = errH1;

  %  % [L2err,H1err] = compute_errors(coordinates,elements3,u);
  % errL2(k) = L2err;
  % errH1(k) = H1err;
end

% compute and print the convergence orders
ordersL2 = log2( L2err(1:end-1) ./ L2err(2:end) );
ordersH1 = log2( H1err(1:end-1) ./ H1err(2:end) );

fprintf('\n   n    L2-error    H1-error   L2-order   H1-order\n');
for j = 1:length(ordersL2)
  fprintf('%4d  %10.3e  %10.3e  %8.3f   %8.3f\n', ...
          ns(j+1), L2err(j+1), H1err(j+1), ordersL2(j), ordersH1(j));
end

% plot on log‑log scale
h = 1./ns;
figure;
loglog(h, L2err, '-o', ...
       h, H1err, '-s', ...
       h, h.^2, '--', ...
       h, h.^1, '-*', ...
       'LineWidth',1.5);
legend('L^2 error','H^1 error','O(h^2)','O(h)','Location','best');
xlabel('h');
ylabel('Error norm');
grid on;

xlabel('h');
ylabel('Error');
title('Error Plot')
set(gca,'XScale','log','YScale','log','LineWidth',1.1);



% 
% %--- Plotting convergence ---
% figure; hold on; box on
% loglog(H,   errL2, 'cy-','LineWidth',1.2);
% loglog(H,   errH1, 'r-','LineWidth',1.2);
% loglog(H,   H.^2,  '--','LineWidth',1.2);
% loglog(H,   H.^1,  '--','LineWidth',1.2);
% legend('L^2','H^1','O(h^2)','O(h)','Location','SouthEast');
% xlabel('h');
% ylabel('Error');
% title('Error Plot')
% set(gca,'XScale','log','YScale','log','LineWidth',1.1);
% 
% %--- Print observed rates ---
% ratesL2 = log(errL2(1:end-1)./errL2(2:end)) ./ log(H(1:end-1)./H(2:end));
% ratesH1 = log(errH1(1:end-1)./errH1(2:end)) ./ log(H(1:end-1)./H(2:end));
% 
% fprintf('\n   h        err_{L2}    rate_{L2}    err_{H1}   rate_{H1}\n');
% for k = 1:Jmax-1
%   fprintf('%6.4f   %8.3e   %7.3f   %8.3e   %7.3f\n', ...
%     H(k), errL2(k), ratesL2(k), errH1(k), ratesH1(k));
% end
