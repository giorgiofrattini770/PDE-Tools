function [A,bv,xnodes,uh] = dtrmixsx(mu, eta, sigma, a, b, delta, beta, fun,N)
% DTRMIXSX: la funzione che risolve il problema ai limiti di
% diffusione-trasporto-reazione con condizioni miste di Dirichlet-Neumann e
% Neumann a sinistra:
%
%         { -mu*u''(x) + eta*u'(x) + sigma(x)*u(x) = f(x) per x in (a,b)
%         { mu*u'(a) = delta
%         { u(b) = beta
%
% dove la derivata prima u'(a) è approssimata con differenze finite
% in avanti (schema accurato di ordine 1).
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
%   beta: il valore della condizione di Dirichlet in u(b)
%   delta: il valore della condizione di Neumann in mu*u(a)
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
%                       uh=(u(0), u(1), u(2), ..., u(N), beta)^T
%
%   Ovvero quindi tale che u(N+1)=u(b)=beta.

h = ( b - a ) / ( N + 1 );
xnodes = linspace( a, b, N + 2 ); 
A = sparse( 2 : N + 1, 2 : N + 1, 2, N + 1, N + 1 ) ...
    + sparse( 2 : N + 1, 1 : N, -1, N + 1, N + 1 ) + sparse( 2 : N, 3 : N + 1, -1, N + 1, N + 1 ); 
A = mu / h^2 * A;
A = A + eta / ( 2 * h ) * ( sparse( 2 : N + 1, 1 : N, -1, N + 1, N + 1 ) ...
                            + sparse( 2 : N, 3 : N + 1, 1, N + 1, N + 1 ) );
A = A + sparse( 2 : N + 1, 2 : N + 1, sigma( xnodes( 2 : end - 1 ) ), N + 1, N + 1 );
A( 1, 1 : 2 ) = mu / h * [ -1 1 ];
bv = ( fun( xnodes( 1 : end - 1 ) ) )';
bv( 1 ) = delta;
bv( end ) = bv( end ) + beta * ( mu / h^2 - eta / ( 2 * h ) );
uh = A \ bv;
uh = [ uh; beta ];
end