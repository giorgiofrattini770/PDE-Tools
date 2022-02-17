function [A,bv,xnodes,uh] = diffreaz(mu, sigma, a, b, alpha, beta, fun,N)
% DIFFREAZ: la funzione che risolve il problema ai limiti di
% diffusione-reazione con condizioni al contorno di Dirichlet, ovvero:
%
%               { -mu*u''(x) + sigma(x)*u(x) = f(x) per x in (a,b)
%               { u(a) = alpha
%               { u(b) = beta
%
% (vedi paragrafo 1 serie 8, in particolare esercizio 1.2)
%
% INPUT:
%   mu: il coefficiente di diffusione per cui è moltiplicato u''(x)
%   sigma: la function handle che esprime il coefficiente di reazione
%   sigma(x) al variare di x. Se tale coefficiente è costante, è
%   sufficiente definire sigma come:
%                       sigma= @(x) costante.*(x==x);
%   a: l'estremo sinistro dell'intervallo di integrazione
%   b: l'estremo destro dell'intervallo di integrazione
%   alpha: il valore della condizione di Dirichlet in u(a)
%   beta: il valore della condizione di Dirichlet in u(b)
%   fun: la forzante definita mediante anonymous function
%   N: il numero di nodi interni (consideriamo N+2 nodi dati da
%   x(j)=x(0)+j*h per j=0,...,N+1 e h=(b-a)/(N+1))
% OUTPUT:
%   A: la matrice CONDENSATA del sistema lineare usato per la risoluzione
%   in formato sparse (in questo caso A è la somma della matrice di
%   diffusione - o  di stiffness - e della matrice diagonale NxN con sigma
%   sulla diagonale)
%   bv: il vettore dei termini noti del sistema lineare CONDENSATO usato
%   per l'approsimazione
%   xnodes: i nodi ottenuti dalla discretizzazione spaziale di passo h
%   uh: il vettore contenente le soluzioni approssimate COMPRESE le
%   condizioni di Dirichlet:
%
%                       uh=(alpha, u(1), u(2), ..., u(N), beta)^T
%
%   Ovvero quindi tale che u(0)=u(x0)=alpha e u(N+1)=u(b)=beta.

h = ( b - a ) / ( N + 1 );
xnodes = linspace( a, b, N + 2 );
A = sparse( 1 : N, 1 : N, 2, N, N ) ...
    + sparse( 2 : N, 1 : N - 1, -1, N, N ) + sparse( 1 : N - 1, 2 : N, -1, N, N );
A = (mu / h^2) * A;
A=A+sparse(diag(sigma(xnodes(2:end-1))));
bv=(fun(xnodes(2:end-1)))';
bv( 1 ) = bv( 1 ) + alpha * mu / h^2;
bv( end ) = bv( end ) + beta * mu / h^2;
uh = A \ bv;
uh = [ alpha; uh; beta ];
end