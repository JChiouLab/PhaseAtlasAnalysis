function varargout = fpeq_model(action, varargin)
%FPEQ_MODEL Full-domain five-species steady-PDE algebra on a unit grid.
% State ordering is species-major: vec(Y.') for Y=[u;q;v;w;p] in R^(5 x nx).

switch action
    case 'laplacian'
        varargout{1} = laplacian(varargin{:});
    case 'weights'
        varargout{1} = trap_weights(varargin{:});
    case 'reaction'
        varargout{1} = reaction(varargin{:});
    case 'reaction_jacobian'
        varargout{1} = reaction_jacobian(varargin{:});
    case 'pde_rhs'
        varargout{1} = pde_rhs(varargin{:});
    case 'pde_jacobian'
        varargout{1} = pde_jacobian(varargin{:});
    case 'augmented'
        [varargout{1:nargout}] = augmented(varargin{:});
    case 'unpack'
        [varargout{1:nargout}] = unpack(varargin{:});
    case 'pack'
        varargout{1} = pack(varargin{:});
    case 'masses'
        varargout{1} = masses(varargin{:});
    otherwise
        error('fpeq_model:unknownAction','Unknown action %s',action);
end
end

function L2 = laplacian(nx)
% Memoized on nx.  pde_rhs calls this once per RHS evaluation, and a time
% integration makes O(10^5) of those; rebuilding the matrix each time meant a
% fresh spdiags plus four sparse element assignments (each a full restructure)
% on every call.  The operator depends only on nx, and no caller mutates the
% returned matrix, so the cached copy is safe to hand out.
persistent NX LL
assert(nx>=5 && nx==round(nx),'nx must be an integer >=5');
if isempty(LL), NX=[]; LL={}; end
hit=find(NX==nx,1);
if ~isempty(hit), L2=LL{hit}; return; end
e=ones(nx,1);
L2=spdiags([e -2*e e],[-1 0 1],nx,nx);
L2(1,1)=-2; L2(1,2)=2;
L2(nx,nx)=-2; L2(nx,nx-1)=2;
L2=L2*(nx-1)^2;
NX(end+1)=nx; LL{end+1}=L2;
if numel(NX)>8, NX=NX(2:end); LL=LL(2:end); end   % bounded: a few meshes per run
end

function w = trap_weights(nx)
dx=1/(nx-1);
w=dx*ones(1,nx);
w([1 end])=dx/2;
end

function R = reaction(Y,th)
u=Y(1,:); q=Y(2,:); v=Y(3,:); w=Y(4,:); p=Y(5,:);
a=th.alpha; b=th.beta; g=th.gamma; d=th.delta; e=th.epsilon;
Ru=(a*u.^2+1).*v-b*u-g*u.*p;
Rq=g*u.*p-d*q;
Rv=-(a*u.^2+1).*v+b*u+e*w;
Rw=d*q-e*w;
Rp=-g*u.*p+d*q;
R=[Ru;Rq;Rv;Rw;Rp];
end

function J = reaction_jacobian(y,th)
u=y(1); v=y(3); p=y(5);
a=th.alpha; b=th.beta; g=th.gamma; d=th.delta; e=th.epsilon;
J=[2*a*u*v-b-g*p, 0, 1+a*u^2, 0, -g*u; ...
   g*p, -d, 0, 0, g*u; ...
   b-2*a*u*v, 0, -(1+a*u^2), e, 0; ...
   0, d, 0, -e, 0; ...
   -g*p, d, 0, 0, -g*u];
end

function rhs = pde_rhs(Y,L,th)
nx=size(Y,2); L2=laplacian(nx);
mu=[th.Dm th.Dm th.Dc th.Dc th.Dc]/(th.f*L^2);
rhs=reaction(Y,th);
for s=1:5
    rhs(s,:)=rhs(s,:)+mu(s)*(L2*Y(s,:).').';
end
end

function A = pde_jacobian(Y,L,th)
nx=size(Y,2); n=5*nx; L2=laplacian(nx);
mu=[th.Dm th.Dm th.Dc th.Dc th.Dc]/(th.f*L^2);
A=sparse(n,n);
for s=1:5
    ii=(s-1)*nx+(1:nx);
    A(ii,ii)=mu(s)*L2;
end
I=zeros(25*nx,1); Jj=I; V=I; c=0;
for k=1:nx
    Jr=reaction_jacobian(Y(:,k),th);
    for s=1:5
        for r=1:5
            c=c+1;
            I(c)=(s-1)*nx+k;
            Jj(c)=(r-1)*nx+k;
            V(c)=Jr(s,r);
        end
    end
end
A=A+sparse(I,Jj,V,n,n);
end

function [F,J,FL,diag] = augmented(z,L,th,nx)
[Y,s1,s2]=unpack(z,nx);
rhs=pde_rhs(Y,L,th);
rhs(1,:)=rhs(1,:)+s1;
rhs(5,:)=rhs(5,:)+s2;
m=masses(Y);
F=[reshape(rhs.',[],1); m(1)-1; m(2)-th.rho];

if nargout>=2
    n=5*nx; A=pde_jacobian(Y,L,th);
    S=sparse(n,2);
    S(1:nx,1)=1;
    S(4*nx+(1:nx),2)=1;
    w=trap_weights(nx);
    C=sparse(2,n);
    C(1,1:nx)=w;
    C(1,nx+(1:nx))=w;
    C(1,2*nx+(1:nx))=w;
    C(1,3*nx+(1:nx))=w;
    C(2,nx+(1:nx))=w;
    C(2,4*nx+(1:nx))=w;
    J=[A S; C sparse(2,2)];
end
if nargout>=3
    L2=laplacian(nx);
    mu=[th.Dm th.Dm th.Dc th.Dc th.Dc]/(th.f*L^2);
    d=zeros(5,nx);
    for s=1:5
        d(s,:)=(-2*mu(s)/L)*(L2*Y(s,:).').';
    end
    FL=[reshape(d.',[],1);0;0];
end
if nargout>=4
    diag=struct('pde_rms',norm(F(1:5*nx))/sqrt(5*nx), ...
        'mass_residual',abs([F(end-1) F(end)]), ...
        'shift_abs',abs([s1 s2]),'minimum',min(Y,[],'all'));
end
end

function z = pack(Y,s1,s2)
if nargin<2, s1=0; end
if nargin<3, s2=0; end
z=[reshape(Y.',[],1);s1;s2];
end

function [Y,s1,s2] = unpack(z,nx)
assert(numel(z)==5*nx+2,'State length does not match nx');
Y=reshape(z(1:5*nx),nx,5).';
s1=z(end-1); s2=z(end);
end

function m = masses(Y)
nx=size(Y,2); w=trap_weights(nx);
m=[w*(sum(Y(1:4,:),1).'), w*((Y(2,:)+Y(5,:)).')];
end

