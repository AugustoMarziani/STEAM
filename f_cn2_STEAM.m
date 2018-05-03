function [result_cn2_teta_wet, Ri, z, T_matrix, Rh_matrix ] = f_cn2_models_matched(matrice_dati)
%RAOB
%La funzione elabora i profili dei vari modelli della Cn2:
%Cn2(z,Theta) in funzione del gradiente della temperatura Theta caso wet
%Funzioni richiamate:
%%AUTORE: Augusto Marziani email: marzianiaugusto@gmail.com



%impostare il vettore in metri dell'altitudine (vettore unidimensionale)





gamma_a=0.0098; % 9.8 °C/km adiabatic lapse rate

% definisco i vettori x,y,z a partire dalla matrice T
[xind, yind, zind] = size(T);

% definisco una matrice Z tridimensionale


Z2 = repmat(z,1,xind,yind);
Z = permute(Z2,[2,3,1]);


% Outer scale Andrew's formula
Lo =5./(1+((Z(:,:,1:end-1)-7500)/2500).^2 );

%imposto i dati meteo:

    %pressione
    P = P;
    
    %temperatura
    T = T;
    
    %valore fittizio per Rh
    rh_test = 55 + (80-55).*rand(1,zind); %genero una colonna variabile di umidità con valore tra 55 e 80
    rh_test = rh_test';
    Rh2 = repmat(rh_test,1,xind,yind);
    Rh = permute(Rh2,[2,3,1]);
        
  
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  
    %Calcolo la temperatura pseudo-potenziale [cherubini et al.]
    H = T + (gamma_a.*Z);
    Hk = H + 273.15 ; % conversione in Kelvin
    
    %Calcolo la temperatura assoluta
    Tk = T + 273.15; % conversione in gradi Kelvin
    
    %vapor saturo
    e_sat = (1.007+3.46*(10^(-6)).*P).*6.1121.*exp(17.502.*T./(240.97+T));
    % T  is the dry bulb temperature expressed in degrees Celsius (°C),
    % P is the absolute pressure expressed in hectopascals,
    %{{e^*}_w}  is the saturated vapor pressure expressed in hectopascals.
    e = (0.01.*Rh) .* e_sat;
    
    %umidità specifica
    Q = e./(1.62.*P);
    
    
    %%Q_mat(:,j) = Q;   ?????
    
    %Temperatura Potenziale Theta [cherubini et al.]
    teta=(Tk).*((1000./P).^0.29); 
    %differenziale temperatura T
    dt_dzeta = diff(Tk,1,3)./diff(Z,1,3);
    %differenziale temperatura potenziale theta
    dteta_dzeta = diff(teta,1,3)./diff(Z,1,3);
    %differenziale temperatura pseudopotenziale H
    dh_dzeta = diff(Hk,1,3)./diff(Z,1,3);
    %differenziale umidità specifica q
    dq_dzeta = diff(Q,1,3)./diff(Z,1,3);
    
    
    %calcolo indice di Richardson
   
%%calcolo di r=mixing ratio of water vapour
r=(0.622.*e)./(P-e);

%%teta_v= teta*(1+0.61r-rl) dove r=mixing ratio of water vapour e rl=mixing
%%ratio of liquid water (rl=0 in aria chiara)
teta_v= teta .*(1+0.61.*r);

%%calcolo della temperatura virtuale
Tv= Tk./ (1-(e./P)*(1-0.622));

U = wind_X; %assumo l'angolo rispetto al nord, il nord è indicato dall'asse Y
V = wind_Y;
dU = (diff(U,1,3)./diff(Z,1,3)).^2;
dV = (diff(V,1,3)./diff(Z,1,3)).^2;
dW = dU + dV;
k2=9.8 ./ teta_v(:,:,2:end);
dteta=(diff(teta_v,1,3)./diff(Z,1,3));
Ri= (k2 .* dteta) ./ dW;
    

    
    %modello Cn2 con temperatura T [cherubini et al.]:
    %caso secco "dry" 
    M_t_dry = (-80*10^(-6).*P(:,:,1:end-1)./(Tk(:,:,1:end-1).^2)) .* dt_dzeta;
    %caso umido "wet"
    M1_t_wet = (1 + (2*4800*1.62.*Q(:,:,1:end-1)./Tk(:,:,1:end-1))).*dt_dzeta;
    M2_t_wet = (4800*1.62).*dq_dzeta;
    M_t_wet = (-80*10^(-6).*P(:,:,1:end-1)./(Tk(:,:,1:end-1).^2)  .* (M1_t_wet - M2_t_wet));
    
    
    %modello Cn2 temperatura pseudo-potenziale H [cherubini et al.]:
    %caso secco "dry"
    M_h_dry = (-80*10^(-6).*P(:,:,1:end-1)./(Tk(:,:,1:end-1).^2)) .* dh_dzeta;
    %caso umido "wet"
    M1_h_wet = (1 + (2*4800*1.62.*Q(:,:,1:end-1)./Tk(:,:,1:end-1))).*dh_dzeta;
    M2_h_wet = (4800*1.62).*dq_dzeta;
    M_h_wet = (-80*10^(-6).*P(:,:,1:end-1)./(Tk(:,:,1:end-1).^2)  .* (M1_h_wet - M2_h_wet));
    
    %modello Cn2 con temperatura potenziale teta [cherubini et al.]:
    %caso secco "dry"    
    M_teta_dry = (-80*10^(-6).*P(:,:,1:end-1)./(Tk(:,:,1:end-1).*teta(:,:,1:end-1))) .* dteta_dzeta;
    %caso umido "wet"
    M1_teta_wet = (1 + (2*4800*1.62.*Q(:,:,1:end-1)./Tk(:,:,1:end-1))).*dteta_dzeta;
    M2_teta_wet = ((4800*1.62)./(Tk(:,:,1:end-1)./teta(:,:,1:end-1))).*dq_dzeta;
    M_teta_wet = (-80*10^(-6).*P(:,:,1:end-1)./(Tk(:,:,1:end-1).^2)  .* (M1_teta_wet - M2_teta_wet));
    
    
    % Cn2_wet(:,j) = (((80*10^-6).*P(2:end))./(T(2:end).^2).*(M_wet).^2; %Formula per il calcolo di Cn2.
    cn2_mod_t_dry = 2.8*1.35.*Lo.^(4/3).*(M_t_dry.^2);
    cn2_mod_t_wet = 2.8*1.35.*Lo.^(4/3).*(M_t_wet.^2);
    
    cn2_mod_h_dry = 2.8*1.35.*Lo.^(4/3).*(M_h_dry.^2);
    cn2_mod_h_wet = 2.8*1.35.*Lo.^(4/3).*(M_h_wet.^2);
    
    cn2_mod_teta_dry = 2.8*1.35.*Lo.^(4/3).*(M_teta_dry.^2);
    cn2_mod_teta_wet = 2.8*1.35.*Lo.^(4/3).*(M_teta_wet.^2);
    
    

  






result_cn2_t_dry = cn2_mod_t_dry;
result_cn2_t_wet = cn2_mod_t_wet;

result_cn2_h_dry = cn2_mod_h_dry;
result_cn2_h_wet = cn2_mod_h_wet;

result_cn2_teta_dry = cn2_mod_teta_dry;
result_cn2_teta_wet = cn2_mod_teta_wet;

Z_cn2 = Z(:,:,1:end-1);

x = result_cn2_teta_wet(1,1,:);
y = Z_cn2(1,1,:);
plot(x(:),y(:));
end
 
