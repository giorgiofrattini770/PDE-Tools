function [A,bv,xnodes,uh,Peh,muh] = difftraspupwind(mu, eta, a, b, alpha, beta, fun,N)
% DIFFTRASP: la funzione che risolve problemi ai limiti di
% diffusione-trasporto con condizioni al contorno di Dirichlet mediante
% tecnica upwind (derivata prima approssimata con differenze finite
% all'indietro o in avanti a seconda del segno di eta):
%
%                       { -mu*u''(x) + eta*u'(x) = f(x) per x in (a,b)
%                       { u(a) = alpha
%                       { u(b) = beta
%
% (vedi paragrafo 2 serie 8, in particolare esercizio 2.1)
% INPUT:
%   mu: il coefficiente di diffusione per cui è moltiplicato u''(x)
%   eta: il coefficiente di trasporto per cui è moltiplicato u'(x)
%   a: l'estremo sinistro dell'intervallo di integrazione
%   b: l'estremo destro dell'intervallo di integrazione
%   alpha: il valore della condizione di Dirichlet in u(a)
%   beta: il valore della condizione di Dirichlet in u(b)
%   fun: la forzante definita mediante anonymous function
%   N: il numero di nodi interni (consideriamo N+2 nodi dati da
%   x(j)=x(0)+j*h per j=0,...,N+1 e h=(b-a)/(N+1))
% OUTPUT:
%   A: la matrice CONDENSATA del sistema lineare usato per la risoluzione
%   in formato sparse
%   bv: il vettore dei termini noti del sistema lineare CONDENSATO usato
%   per l'approsimazione
%   xnodes: i nodi ottenuti dalla discretizzazione spaziale di passo h
%   uh: il vettore contenente le soluzioni approssimate COMPRESE le
%   condizioni di Dirichlet:
%                       uh=(alpha, u(1), u(2), ..., u(N), beta)^T
%   Ovvero quindi tale che u(0)=u(x0)=alpha e u(N+1)=u(b)=beta.
%   Peh: numero di Peclet locale associato al problema di
%   diffusione-trasporto
%   muh: diffusività artificiale associata al problema numerico con tecnica
%   upwind

h = ( b - a ) / ( N + 1 );
xnodes = linspace( a, b, N + 2 );
Peh = abs( eta ) * h / ( 2 * mu );
muh = mu * ( 1 + Peh );
A = sparse( 1 : N, 1 : N, 2, N, N ) ...
    + sparse( 2 : N, 1 : N - 1, -1, N, N ) + sparse( 1 : N - 1, 2 : N, -1, N, N ); 
A = muh / h^2 * A;
A = A + eta / ( 2 * h ) * ( sparse( 2 : N, 1 : N - 1, -1, N, N ) ...
                            + sparse( 1 : N - 1, 2 : N, 1, N, N ) );
bv = ( fun( xnodes( 2 : end - 1 ) ) )';
bv( 1 ) = bv( 1 ) + alpha * ( muh / h^2 + eta / ( 2 * h ) );
bv( end ) = bv( end ) + beta * ( muh / h^2 - eta / ( 2 * h ) );
uh = A \ bv;
uh = [ alpha; uh; beta ];
end