% RiskBrowser.m
disp('********************************************************************')
disp('*   Exact Risk Analysis in Wavelet Regression                      *');
disp('********************************************************************')
risk_global;
%RiskBrowserIntro

back_color = [0.9 0.9 0.9];
global back_color
num_menus = 5;
  menu_names = [ ...
    'Data______'; ...
    'Signals___'; ...
    'Settings__'; ...
    'Actions___'; ...
    'Plots_____'; ...
               ];

% To add menu items, increase the appropriate number in the array
% below, and add the corresponding 'entries' and 'call backs':

num_menu_items = [ 1 10 0 5 0];

% -- Data for Data menu # 1 --- %

Data_______entries = [ ...
         'User Defined' ; ...
                     ];

Data_______callbacks = [ ...
         'Gen_data( 0 );' ; ...
                       ];

% -- Data for Signals menu #2 --- %

Signals____entries = [ ...
         'Step             ' ; ... 
         'Wave             ' ; ...
         'Blip             ' ; ...
         'Blocks           ' ; ...
         'Bumps            ' ; ...
         'HeaviSine        ' ; ...
         'Doppler          ' ; ...
         'Angles           ' ; ...
         'Parabolas        ' ; ...
         'TSh Sine         ' ; ...
                   ];


% -- Data for Settings menu #3 --- %

Settings___entries = [ ...
           'Signal Length     ' ; ...
           'Noise Type        ' ; ...
           'Noise Level       ' ; ...
           'Wavelet           ' ; ...
           'Threshold selector' ; ...
           'Plot options      ' ; ...
                    ];

% -- Data for Actions menu #3 --- %

Actions____entries = [ ...
           'Add Noise             ' ; ...
           'Fourier Transform     ' ; ...
           'Wavelet Transform     ' ; ...
           'DeNoise/Compress      ' ; ...
           'Quit                  ' ; ...
                     ];


Plots______entries = [ ...
         'Wavelet and Fourier Transforms' ; ...
         'Power Remaining               ' ; ...
         'Exact Risk Plots              ' ; ...
                     ];

%      Put in the other submenus now
%      Settings menu


        %   set up default values
%
        n = 256;
        x_length = n;
        x_name = ' ';
        if isempty(Wav_type),
          Wav_type = 'Haar';
        end
        if isempty(threshtype),
          threshtype = 'Hard';
        end
        if isempty(noiseamp),
         noiseamp = 0.02;       % Low Noise
        end
        if isempty(noisetype),
          noisetype = 'Normal';
        end;
        Empty_Data = 1;
        signal_name = '';
        panel = n/2;
        ylim  = [-1 1];
        x_use  = zeros(1,n);
        sigchoice = zeros(1,10);
%
 
    
    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:42 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
