function [A,bv,xnodes,uh] = dtrmixdxord2(mu, eta, sigma, a, b, alpha, gamma, fun,N)
% DTRMIXDXORD2: la funzione che risolve il problema ai limiti di
% diffusione-trasporto-reazione con condizioni miste di Dirichlet-Neumann e
% Neumann a destra:
%
%         { -mu*u''(x) + eta*u'(x) + sigma(x)*u(x) = f(x) per x in (a,b)
%         { u(a) = alpha
%         { mu*u'(b) = gamma
%
% dove la derivata prima u'(b) è approssimata con uno schema decentrato per
% la derivata prima di ordine 2.
% (vedi paragrafo 3 serie 8, in particolare esercizio 3.1)
%
% INPUT:
%   mu: il coefficiente di diffusione per cui è moltiplicato u''(x)
%   eta: il coefficiente di trasporto per cui è moltiplicato u'(x)
%   sigma: la function handle che esprime il coefficiente di reazione
%   sigma(x) al variare di x. Se tale coefficiente è costante, è
%   sufficiente definire sigma come:
%                       sigma= @(x) costante.*(x==x);
%   a: l'estremo sinistro dell'intervallo di integrazione
%   b: l'estremo destro dell'intervallo di integrazione
%   alpha: il valore della condizione di Dirichlet in u(a)
%   gamma: il valore della condizione di Neumann in mu*u(b)
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
%                       uh=(alpha, u(1), u(2), ..., u(N), u(N+1))^T
%
%   Ovvero quindi tale che u(0)=u(x0)=alpha

h = ( b - a ) / ( N + 1 );
xnodes = linspace( a, b, N + 2 );
A = sparse( 1 : N, 1 : N, 2, N + 1, N + 1 ) ...
    + sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) + sparse( 1 : N, 2 : N + 1, -1, N + 1, N + 1 ); 
A = mu / h^2 * A;
A = A + eta / ( 2 * h ) * ( sparse( 2 : N, 1 : N - 1, -1, N + 1, N + 1 ) ...
                            + sparse( 1 : N, 2 : N + 1, 1, N + 1, N + 1 ) );
A = A + sparse( 1 : N, 1 : N, sigma( xnodes( 2 : end - 1 ) ), N + 1, N + 1 );
A( end, N - 1 : N + 1 ) = mu / ( 2 * h ) * [ 1 -4 3 ];
bv = ( fun( xnodes( 2 : end ) ) )';
bv( 1 ) = bv( 1 ) + alpha * ( mu / h^2 + eta / ( 2 * h ) );
bv( end ) = gamma;
uh = A \ bv;
uh = [ alpha; uh ];
end