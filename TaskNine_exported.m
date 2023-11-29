classdef TaskNine_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ImageProcessing                matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        ImageViewerTab                 matlab.ui.container.Tab
        BrowseButton                   matlab.ui.control.Button
        PathEditField                  matlab.ui.control.EditField
        WidthLabel                     matlab.ui.control.Label
        HeightLabel                    matlab.ui.control.Label
        SizeinbitsLabel                matlab.ui.control.Label
        BitDepthLabel                  matlab.ui.control.Label
        ImageColorLabel                matlab.ui.control.Label
        WidthEditField                 matlab.ui.control.EditField
        HeightEditField                matlab.ui.control.EditField
        SizeInBitsEditField            matlab.ui.control.EditField
        BitDepthEditField              matlab.ui.control.EditField
        ImageColorEditField            matlab.ui.control.EditField
        ModalityLabel                  matlab.ui.control.Label
        PatientNameLabel               matlab.ui.control.Label
        PatientAgeLabel                matlab.ui.control.Label
        BodyPartExaminedLabel          matlab.ui.control.Label
        ModalityEditField              matlab.ui.control.EditField
        PatientNameEditField           matlab.ui.control.EditField
        PatientAgeEditField            matlab.ui.control.EditField
        BodyPartExaminedEditField      matlab.ui.control.EditField
        ImageDisplay                   matlab.ui.control.UIAxes
        InterpolationTab               matlab.ui.container.Tab
        InterpolateButton              matlab.ui.control.Button
        ZoomingFactorLabel             matlab.ui.control.Label
        ZoomingFactorEditField         matlab.ui.control.NumericEditField
        InterpolationTypeDropDown      matlab.ui.control.DropDown
        InterpolationPanel             matlab.ui.container.Panel
        InterpolatedImage              matlab.ui.control.UIAxes
        OriginalPanel                  matlab.ui.container.Panel
        OriginalImage                  matlab.ui.control.UIAxes
        OriginalImageinGrayscaleLabel  matlab.ui.control.Label
        InterpolatedImageLabel         matlab.ui.control.Label
        OriginalWidthLabel             matlab.ui.control.Label
        OriginalHeightLabel            matlab.ui.control.Label
        NewWidthLabel                  matlab.ui.control.Label
        NewHeightLabel                 matlab.ui.control.Label
        OriginalWidthEditField         matlab.ui.control.EditField
        OriginalHeightEditField        matlab.ui.control.EditField
        NewWidthEditField              matlab.ui.control.EditField
        NewHeightEditField             matlab.ui.control.EditField
        info                           matlab.ui.control.Label
        AffineTransformationTab        matlab.ui.container.Tab
        ConstructImageButton           matlab.ui.control.Button
        AngleSlider                    matlab.ui.control.Slider
        AngleSliderLabel               matlab.ui.control.Label
        InterpolationTypeOfRotatedImageDropDown  matlab.ui.control.DropDown
        InterpolationTypeLabel         matlab.ui.control.Label
        ShearingAngleLabel             matlab.ui.control.Label
        Shear45Button                  matlab.ui.control.Button
        ShearNegative45Button          matlab.ui.control.Button
        AngleLabel                     matlab.ui.control.Label
        DirectionLabel                 matlab.ui.control.Label
        AngleEditField                 matlab.ui.control.EditField
        DirectionEditField             matlab.ui.control.EditField
        RotateWidthLabel               matlab.ui.control.Label
        RotateHeightLabel              matlab.ui.control.Label
        RotateWidthEditField           matlab.ui.control.EditField
        RotateHeightEditField          matlab.ui.control.EditField
        TImage                         matlab.ui.control.UIAxes
        HistogramTab                   matlab.ui.container.Tab
        HistogramBrowseButton          matlab.ui.control.Button
        EqualizedImageLabel            matlab.ui.control.Label
        OriginalImageLabel             matlab.ui.control.Label
        UnequalizedDisplay             matlab.ui.control.UIAxes
        EqualizedDisplay               matlab.ui.control.UIAxes
        UnequalizedHistogram           matlab.ui.control.UIAxes
        EqualizedHistogram             matlab.ui.control.UIAxes
        SpatialFilteringTab            matlab.ui.container.Tab
        SpatialFilterBrowseButton      matlab.ui.control.Button
        BlurredNoisyImageLabel         matlab.ui.control.Label
        KernelsizeLabel                matlab.ui.control.Label
        FactorKLabel                   matlab.ui.control.Label
        KernelSizeSpinner              matlab.ui.control.Spinner
        EnhanceButton                  matlab.ui.control.Button
        RemoveNoiseButton              matlab.ui.control.Button
        EnhancedImageLabel             matlab.ui.control.Label
        FactorKEditField               matlab.ui.control.NumericEditField
        UnfilteredImageDisplay         matlab.ui.control.UIAxes
        Blurred_NoisyDisplay           matlab.ui.control.UIAxes
        EnhancedImageDisplay           matlab.ui.control.UIAxes
        FourierITab                    matlab.ui.container.Tab
        FourierBrowseButton            matlab.ui.control.Button
        LogDropDown                    matlab.ui.control.DropDown
        MagnitudeLabel                 matlab.ui.control.Label
        PhaseLabel                     matlab.ui.control.Label
        OriginalImageLabel_2           matlab.ui.control.Label
        OriginalImageDisplay           matlab.ui.control.UIAxes
        MagnitudeDisplay               matlab.ui.control.UIAxes
        PhaseDisplay                   matlab.ui.control.UIAxes
        FourierIITab                   matlab.ui.container.Tab
        FourierIIBrowseButton          matlab.ui.control.Button
        FourierKernelsizeLabel         matlab.ui.control.Label
        FourierKernelSizeSpinner       matlab.ui.control.Spinner
        FilterButton                   matlab.ui.control.Button
        FilteredImageLabel             matlab.ui.control.Label
        DifferenceLabel                matlab.ui.control.Label
        ImageDisplaay                  matlab.ui.control.UIAxes
        FourierFilterDisplay           matlab.ui.control.UIAxes
        DifferenceDisplay              matlab.ui.control.UIAxes
        RemovePatternTab               matlab.ui.container.Tab
        RemovePatternButton            matlab.ui.control.Button
        PatternedImageLabel            matlab.ui.control.Label
        MaskingLabel                   matlab.ui.control.Label
        PatternRemovedLabel            matlab.ui.control.Label
        Mask                           matlab.ui.control.UIAxes
        ImageWithoutPattern            matlab.ui.control.UIAxes
        ImageWithPattern               matlab.ui.control.UIAxes
        NoiseTab                       matlab.ui.container.Tab
        CreateImageButton              matlab.ui.control.Button
        NoiseTypeDropDown              matlab.ui.control.DropDown
        NoiseTypeLabel                 matlab.ui.control.Label
        ConfirmRegionofInterestButton  matlab.ui.control.Button
        MeanLabel                      matlab.ui.control.Label
        MeanEditField                  matlab.ui.control.EditField
        VarianceLabel                  matlab.ui.control.Label
        VarianceEditField              matlab.ui.control.EditField
        CreatedImage                   matlab.ui.control.UIAxes
        ROIHistogram                   matlab.ui.control.UIAxes
        BackProjectionTab              matlab.ui.container.Tab
        StartButton                    matlab.ui.control.Button
        ScheppLoganPhantomLabel        matlab.ui.control.Label
        SinogramLabel                  matlab.ui.control.Label
        LaminogramLabel                matlab.ui.control.Label
        LaminogramDropDown             matlab.ui.control.DropDown
        ManualButton                   matlab.ui.control.Button
        BuiltinButton                  matlab.ui.control.Button
        ScheppLoganPhantom             matlab.ui.control.UIAxes
        Sinogram                       matlab.ui.control.UIAxes
        Laminogram                     matlab.ui.control.UIAxes
        ImageMorphologyTab             matlab.ui.container.Tab
        DisplayButton                  matlab.ui.control.Button
        FingerprintImageLabel          matlab.ui.control.Label
        MorphologicalProcessingLabel   matlab.ui.control.Label
        DilationButton                 matlab.ui.control.Button
        ErosionButton                  matlab.ui.control.Button
        ClosingButton                  matlab.ui.control.Button
        OpeningButton                  matlab.ui.control.Button
        Button                         matlab.ui.control.Button
        Button_2                       matlab.ui.control.Button
        Button_3                       matlab.ui.control.Button
        Button_4                       matlab.ui.control.Button
        RotateWidthLabel_2             matlab.ui.control.Label
        RotateHeightLabel_2            matlab.ui.control.Label
        FingerprintImage               matlab.ui.control.UIAxes
        MorphologicalImage             matlab.ui.control.UIAxes
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: BrowseButton
        function BrowseButtonPushed(app, event)
 %Open dialog box for browsing
 global Filename;
 global Pathname;
 global indx;
            Filename=0;
            [Filename, Pathname,indx] = uigetfile({'*.jpg';'*.jpeg';'*.bmp';'*.dcm';'*.tif';'*.png'}, 'All Files');

            if Filename==0
                return;
            end
            %if extensions jpg or bmp
             if (indx==1)||(indx==2)||(indx==3)||(indx==5)||(indx==6)
                fullname=[Pathname,Filename];
                app.PathEditField.Value=fullname;
                app.ModalityLabel.Visible='off';
                app.ModalityEditField.Visible='off';
                app.PatientNameLabel.Visible='off';
                app.PatientNameEditField.Visible='off';
                app.PatientAgeLabel.Visible='off';
                app.PatientAgeEditField.Visible='off';
                app.BodyPartExaminedLabel.Visible='off';
                app.BodyPartExaminedEditField.Visible='off';
                
                %try & catch error due to corrupted image
                %Display image with its info
                try
                    ImageFile=imread(fullname);
                    imshow(ImageFile,'Parent' ,app.ImageDisplay);
                    [row,col,channels]=size(ImageFile);
                    ImageData=imfinfo(fullname);
                    app.WidthEditField.Value=num2str(col);
                    app.HeightEditField.Value=num2str(row);
                    app.BitDepthEditField.Value=num2str(ImageData.BitDepth);
                    SizeinBits=col*row*ImageData.BitDepth;
%                   
                    
                    
                    app.SizeInBitsEditField.Value=num2str(SizeinBits);
                    if (ImageData.BitDepth==1 && channels==1 )
                        app.ImageColorEditField.Value='Binary';
                    elseif (ImageData.BitDepth==8 && channels==1)
                        app.ImageColorEditField.Value='Grayscale';
                    else
                       app.ImageColorEditField.Value='RGB'; 
                    end
       
                catch 
                    warndlg('Corrupted Image! Please rebrowse.','Warning');
                end
              
             elseif (indx==4)
                fullname=[Pathname,Filename];
                app.PathEditField.Value=fullname;
                app.ModalityLabel.Visible='on';
                app.ModalityEditField.Visible='on';
                app.PatientNameLabel.Visible='on';
                app.PatientNameEditField.Visible='on';
                app.PatientAgeLabel.Visible='on';
                app.PatientAgeEditField.Visible='on';
                app.BodyPartExaminedLabel.Visible='on';
                app.BodyPartExaminedEditField.Visible='on';
                
                %try & catch error due to corrupted dicom image
                %Display dicom image with its image and header info
                try
                    ImageFile=dicomread(fullname);
                    app.ImageDisplay.Visible='on';
                    imshow(ImageFile,'Parent',app.ImageDisplay);
                    DicomImgData=dicominfo(fullname);
                    app.WidthEditField.Value=num2str(DicomImgData.Width);
                    app.HeightEditField.Value=num2str(DicomImgData.Height);
                    app.BitDepthEditField.Value=num2str(DicomImgData.BitDepth);
                    SizeinBitsdcm=double(DicomImgData.Width)*double(DicomImgData.Height)*double(DicomImgData.BitDepth);
                    app.SizeInBitsEditField.Value=num2str(SizeinBitsdcm);
                    app.ImageColorEditField.Value=(DicomImgData.ColorType);
                    app.PathEditField.Value=fullname;
                    app.ModalityEditField.Value=char(DicomImgData.Modality);

                    try
                        app.PatientAgeEditField.Value=num2str(DicomImgData.PatientAge); 
                    catch 
                        app.PatientAgeEditField.Value='Not Found';
                    end
                    
                    
                    try
                        app.BodyPartExaminedEditField.Value=char(DicomImgData.BodyPartExamined);
                    catch 
                        app.BodyPartExaminedEditField.Value='Not Found';
                    end
                    
 
                    try
                        app.PatientNameEditField.Value=struct2array(DicomImgData.PatientName);
                    catch 
                        app.PatientNameEditField.Value='Not Found';   
                    end
   
                catch 
                    Filename=0;
                   warndlg('Corrupted Dicom Image! Please rebrowse.','Warning');
                end
          
                
             else
                 Filename=0;
                warndlg('Please choose an image file');
             end
        end

        % Button pushed function: InterpolateButton
        function InterpolateButtonPushed(app, event)
global Filename;
global Pathname;
global indx;
global ImageFile;
            
            %Checks if there is an image chosen
            try
                fullname=[Pathname,Filename];
                % converts dicom to jpg
                if (indx==4)
                    dcmfile = dicomread(fullname);
                    dcmImage = uint8(255 * mat2gray(dcmfile));
                    folder='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task2-MariamOkasha';
                    imwrite(dcmImage,fullfile(folder,'dcmImage.jpg'));
                    fullname='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task2-MariamOkasha/dcmImage.jpg';
                end
                
                ImageFile=imread(fullname); 
                [row,col,channels]=size(ImageFile);
                %displays original dimensions
                app.OriginalWidthEditField.Value=num2str(col);
                app.OriginalHeightEditField.Value=num2str(row);
                        
                %Converts RGB to greyscale
                if channels==3
                    ImageFile= rgb2gray(ImageFile); 
                end
                
                app.OriginalImage.Position=[0,0,row,col];
                imshow(ImageFile,[],'Parent' ,app.OriginalImage); 
                
                app.ZoomingFactorLabel.Visible='on';
                app.ZoomingFactorEditField.Visible='on';
                app.InterpolationTypeDropDown.Visible='on';

            catch 
                warndlg('Please choose an image from "Image Viewer" Tab.','Warning');
            end       
            Filename=0;
        end

        % Value changed function: InterpolationTypeDropDown
        function InterpolationTypeDropDownValueChanged(app, event)
global ImageFile;        
            InterpolationType = app.InterpolationTypeDropDown.Value;
            
            [row, col] = size(ImageFile);

            ZoomingFactor=app.ZoomingFactorEditField.Value; 
            
            
            NewRowCount=round(ZoomingFactor*row);
            NewColCount=round(ZoomingFactor*col);
            
            app.NewWidthEditField.Value=num2str(NewColCount);
            app.NewHeightEditField.Value=num2str(NewRowCount);
            
            if (strcmp(InterpolationType,'Nearest-Neighbour Interpolation')==1)
                %Loop of row
                RowRange = 1 : NewRowCount;
                AlterRow = round(RowRange * (row - 1) / (NewRowCount - 1) + (NewRowCount - row) / (NewRowCount- 1)-0.1 );
                
                %Loop of col
                ColRange = 1 : NewColCount;
                AlterCol = round(ColRange * (col - 1) / (NewColCount - 1) + (NewColCount - col) / (NewColCount - 1)-0.1 ); 
                
                nearestinterpol(RowRange,ColRange) = ImageFile(AlterRow,AlterCol);
                app.InterpolatedImage.Position=[0,0,NewRowCount,NewColCount];
                imshow(nearestinterpol,[],'Parent' ,app.InterpolatedImage);
                app.InterpolatedImageLabel.Text = 'Nearest-Neighbour Interpolation';
                
                
                
            elseif (strcmp(InterpolationType,'Linear Interpolation')==1) 
               
                
                
                LinInterpol = zeros(NewRowCount,NewColCount);
              
                for i=1:NewRowCount
                    %row index of output 
                     y = (NewRowCount-i)/(NewRowCount-1)+row*((i-1)/(NewRowCount-1));
                   
                     for j=1:NewColCount 
                        %col index of output
                        x = (NewColCount-j)/(NewColCount-1)+col*((j-1)/(NewColCount-1));
                 
                        %Any values out of acceptable range
                        
                        x(x > col - 0.00001) = col - 0.00001;
                        x1 = floor(x);
                        x2 = x1 + 1;
                        
                        y(y > row - 0.00001) = row - 0.00001; 
                        y1 = floor(y);
                        y2 = y1 + 1;
                        
                        %// 4 Neighboring Pixels
                        NP1 = ImageFile(y1,x1);
                        NP2 = ImageFile(y1,x2);
                        NP3 = ImageFile(y2,x1); 
                        NP4 = ImageFile(y2,x2);
                        %// 4 Pixels Weights
                        PW1 = (y2-y)*(x2-x);
                        PW2 = (y2-y)*(x-x1);
                        PW3 = (x2-x)*(y-y1);
                        PW4 = (y-y1)*(x-x1);
                        LinInterpol(i,j) = PW1 * NP1 + PW2 * NP2 + PW3 * NP3 + PW4 * NP4;
                    end
                end
                app.InterpolatedImage.Position=[0,0,NewRowCount,NewColCount];
                
                
                
                imshow(LinInterpol,[],'Parent' ,app.InterpolatedImage);
                app.InterpolatedImageLabel.Text = 'Linear Interpolation';
                
            end

            
        end

        % Button pushed function: ConstructImageButton
        function ConstructImageButtonPushed(app, event)
            global TShape;
            app.AngleSlider.Value=0;
            app.InterpolationTypeOfRotatedImageDropDown.Value = app.InterpolationTypeOfRotatedImageDropDown.Items(1);
            
            TShape = zeros(128,128, 'uint8');
            TShape(29:49, 29:99) = 255;
            TShape(49:99, 54:74) = 255;
           
            imshow(TShape,'Parent' ,app.TImage); 
        end

        % Value changed function: AngleSlider
        function AngleSliderValueChanged(app, event)
            global TShape;
            global RotatedTShape;
            global angle;
            global xo;
            global yo;
            global midx;
            global midy;
            
            app.InterpolationTypeOfRotatedImageDropDown.Value=app.InterpolationTypeOfRotatedImageDropDown.Items(1);
            angle = app.AngleSlider.Value;
            
            [row,col]= size(TShape);

           
            NewRowCount=ceil(row*abs(cosd(angle))+col*abs(sind(angle)));                     
            NewColCount=ceil(row*abs(sind(angle))+col*abs(cosd(angle))); 
            
            app.RotateHeightEditField.Value=num2str(NewRowCount);
            app.RotateWidthEditField.Value=num2str(NewColCount);
            
            RotatedTShape=uint8(zeros([NewRowCount NewColCount]));
            
            %calculating center of original a image
            xo=ceil(row/2);                                                            
            yo=ceil(col/2);
            
            %calculating center of final image
            midx=ceil((size(RotatedTShape,1))/2); 
            midy=ceil((size(RotatedTShape,2))/2); 
            
            
            m=1;
            for i=1:size(TShape,1)
                for j=1:size(TShape,2)
                    %convert from cartesian to polar
                    [theta,r]=cart2pol(i-xo,j-yo);
                    %Convert from radians to degree and add the degree value
                    theta1=radtodeg(theta)+angle;
                    %Convert from degree to radians
                    theta=degtorad(theta1);
                    %Convert to Cartesian Co-ordinates
                    [x,y]=pol2cart(theta,r);
                    x1(m)=round(x+midx);
                    x2(m)=round(y+midy);
                  
                    m=m+1;
                end
            end
            %check whether the values are within the image size.
            x1(x1 < 1)=1;
            x2(x2 < 1)=1;
            
            n=1;
            for i=1:size(TShape,1)
                for j=1:size(TShape,2)
                    RotatedTShape(x1(n),x2(n),:)=TShape(i,j,:);
                   
                    n=n+1;
                end 
            end 
            
            imshow(RotatedTShape,'Parent' ,app.TImage);
            app.AngleEditField.Value=num2str(abs(angle));
            if (angle>0)
                app.DirectionEditField.Value="Anti-Clockwise";
            elseif (angle<0)
                app.DirectionEditField.Value="Clockwise";
            else
                app.DirectionEditField.Value="Zero";
            end
           
        end

        % Value changed function: 
        % InterpolationTypeOfRotatedImageDropDown
        function InterpolationTypeOfRotatedImageDropDownValueChanged(app, event)
            
            global TShape;  
            global angle;


            
            InterpolationTypeOfRotatedImage = app.InterpolationTypeOfRotatedImageDropDown.Value;
            [row, col] = size(TShape); 
            xo=ceil(row/2);                                                            
            yo=ceil(col/2);
            if (strcmp(InterpolationTypeOfRotatedImage,'Select Interpolation Method')==1)
              app.AngleSliderValueChanged();
              
            elseif (strcmp(InterpolationTypeOfRotatedImage,'Nearest-Neighbour Interpolation')==1)
                if (angle==0)
                   imshow(TShape,'Parent' ,app.TImage);
                else
                  
                   NewRowCountt=ceil(row*abs(cosd(angle))+col*abs(sind(angle))); 
                   NewColCountt=ceil(row*abs(sind(angle))+col*abs(cosd(angle))); 
                    
                   app.RotateHeightEditField.Value=num2str(NewRowCountt);
                   app.RotateWidthEditField.Value=num2str(NewColCountt);                
                    
                    
                   
                    RotatedTShapee = uint8(zeros(NewRowCountt, NewColCountt));
                    for i=1:NewRowCountt
                        for j=1:NewColCountt                                                    
                    
                             x= (i-NewRowCountt/2)*cosd(angle) + (j-NewColCountt/2)*sind(angle);              
                             y= -(i-NewRowCountt/2)*sind(angle) + (j-NewColCountt/2)*cosd(angle);                           
                             x=round(x+xo); 
                             y=round(y+yo); 
                    
                             if (x>=1 && y>=1 && x<=size(TShape,1) &&  y<=size(TShape,2) ) 
                                  RotatedTShapee(i,j,:)=TShape(x,y,:); 
                             end
                        end
                    end
                     imshow(RotatedTShapee,'Parent' ,app.TImage);
                end
               
                
                
            elseif (strcmp(InterpolationTypeOfRotatedImage,'Linear Interpolation')==1)
                if (angle==0)
                    imshow(TShape,'Parent' ,app.TImage);
                else
                  
                   NewRowCount=ceil(row*abs(cosd(angle))+col*abs(sind(angle)));          
                   NewColCount=ceil(row*abs(sind(angle))+col*abs(cosd(angle)));           
                    
                   app.RotateHeightEditField.Value=num2str(NewRowCount);
                   app.RotateWidthEditField.Value=num2str(NewColCount);                
                    
                    
                    % Original needs to be the image I
                    RotatedTShapeee = uint8(zeros(NewRowCount, NewColCount));
                   
                    
                    for i = 1:NewRowCount
                        for j = 1:NewColCount
                             yo = (j - NewColCount/2) * cosd(angle) - (i - NewRowCount/2)*sind(angle)+64 ; 
                             xo = (j - NewColCount/2) * sind(angle) + (i - NewRowCount/2)*cosd(angle) +64; 
                    
                             x1 = floor(xo);  
                             y1 = floor(yo); 
                             a = xo - x1; 
                             b = yo - y1;
                             
                             % check if the coordinate is in the original image to prevent errors
                             if x1 > 0  && y1 > 0 && y1 <= col && x1 <= row
                                 RotatedTShapeee(i, j) = RotatedTShapeee(i, j) + (1-a)*(1-b)*TShape(x1,y1);
                             end
                             if x1 > 0  && y1+1 > 0 && y1+1 <= col && x1 <= row
                                 RotatedTShapeee(i, j) = RotatedTShapeee(i, j) + (1-a)*b*TShape(x1,y1+1);
                             end
                             if x1+1 > 0  && y1 > 0 && y1 <= col && x1+1 <= row
                                 RotatedTShapeee(i, j) = RotatedTShapeee(i, j) + a*(1-b)*TShape(x1+1,y1);
                             end
                             if x1+1 > 0  && y1+1 > 0 && y1+1 <= col && x1+1 <= row
                                 RotatedTShapeee(i, j) = RotatedTShapeee(i, j) + a*b*TShape(x1+1,y1+1);
                             end
                        end
                    end
                    imshow(RotatedTShapeee,'Parent' ,app.TImage);
                end
           
            end        
        end

        % Button pushed function: Shear45Button
        function Shear45ButtonPushed(app, event)
              global TShape;
              [row,col]= size(TShape);
              newCol=col+row-1;
              
             shear=zeros([row newCol]);
             
              for i=1:row
                 shearedCol=newCol-i+1;
                 
                 for j=shearedCol-col+1:shearedCol
                    
                     shear(i,j)=TShape(i,(j-col+i));
                    
                 end
             end
             imshow(shear,'Parent' ,app.TImage);
        end

        % Button pushed function: ShearNegative45Button
        function ShearNegative45ButtonPushed(app, event)

              global TShape;
              [row,col]= size(TShape);
              newCol=col+row-1;
              
             shear=zeros([row newCol]);
             for i=1:row
                 shearedCol=i+col-1;
                 for j=i:shearedCol
                    
                     shear(i,j)=TShape(i,(j-i+1));
                    
                 end
             end
             imshow(shear,'Parent' ,app.TImage);
        end

        % Button pushed function: HistogramBrowseButton
        function HistogramBrowseButtonPushed(app, event)
            [Filename, Pathname,indx] = uigetfile({'*.jpg';'*.jpeg';'*.bmp';'*.dcm';'*.tif';'*.png'}, 'All Files');
             
            if Filename==0
                return;    
            end
            
            fullname=[Pathname,Filename];
            if (indx==4)
                dcmfile = dicomread(fullname);
                dcmImage = uint8(255 * mat2gray(dcmfile));
                folder='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task4-MariamOkasha/Images';
                imwrite(dcmImage,fullfile(folder,'dcmImage.jpg'));
                fullname='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task4-MariamOkasha/Images/dcmImage.jpg';
            end
       
            try
                ImageFile=imread(fullname);
                [row, col, channels] = size(ImageFile);
                  
                if channels==3
                    ImageFile= rgb2gray(ImageFile);
                end
                imshow(ImageFile,'Parent' ,app.UnequalizedDisplay);
                
                %%%%%%%%%%%%%%%%%%%%%%normalized histogram%%%%%%%%%%%%%%%%%%%%%%%
                
                freq=zeros(256,1);
                normalized=zeros(256,1);
           
                for i=1:row
                    for j=1:col 
    	                pix=ImageFile(i,j);
    	                freq(pix+1)=freq(pix+1)+1;                               
                        normalized(pix+1)=freq(pix+1)/(row*col);                        
                    end
                end
                
                bar(app.UnequalizedHistogram,normalized);
                
                %%%%%%%%%%%%%%%%%%%%%%equalization%%%%%%%%%%%%%%%%%%%%%%%
                
                EqualizedImageFile = uint8(zeros(row,col));
                cdf = zeros(256,1);
                cumulative = zeros(256,1);
                output = zeros(256,1);
                
                s =0 ;
                intensityLevel = 255;
                
                for i=1:size(normalized)
                   s=s+freq(i);                
                   cumulative(i)=s;                
                   cdf(i)=cumulative(i)/(row*col);                
                   output(i)=round(cdf(i)*intensityLevel);                
                end
                
                for i=1:row               
                    for j=1:col
                         EqualizedImageFile(i,j)=output(ImageFile(i,j)+1);                
                    end                
                end  
                
                imshow(EqualizedImageFile,'Parent' ,app.EqualizedDisplay);
                
                %%%%%%%%%%%%%%%%%%%%%%equalization histogram%%%%%%%%%%%%%%%%%%%%%%%
                
                [eqrow, eqcol] = size(ImageFile);
                eqfreq=zeros(256,1);
                eqnormalized=zeros(256,1);
           
                for i=1:eqrow
                    for j=1:eqcol
    	                eqpix=EqualizedImageFile(i,j);
    	                eqfreq(eqpix+1)=eqfreq(eqpix+1)+1;                               
                        eqnormalized(eqpix+1)=eqfreq(eqpix+1)/(eqrow*eqcol);                        
                    end
                end
                
                bar(app.EqualizedHistogram,eqnormalized); 
                
                
            catch 
                warndlg('Corrupted Image! Please rebrowse.','Warning');
            end
        end

        % Button pushed function: SpatialFilterBrowseButton
        function SpatialFilterBrowseButtonPushed(app, event)
          global ImageFile;
            [Filename, Pathname,indx] = uigetfile({'*.jpg';'*.jpeg';'*.bmp';'*.dcm';'*.tif';'*.png'}, 'All Files');
             
            if Filename==0
                return;    
            end
            
            fullname=[Pathname,Filename];
            if (indx==4)
                dcmfile = dicomread(fullname);
                dcmImage = uint8(255 * mat2gray(dcmfile));
                folder='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task4-MariamOkasha/Images';
                imwrite(dcmImage,fullfile(folder,'dcmImage.jpg'));
                fullname='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task4-MariamOkasha/Images/dcmImage.jpg';
            end
       
            try
                ImageFile=imread(fullname);
                [row, col, channels] = size(ImageFile);
                  
                if channels==3
                    ImageFile= rgb2gray(ImageFile);
                end
             imshow(ImageFile,'Parent' ,app.UnfilteredImageDisplay);
             app.RemoveNoiseButton.Visible = 'on';
             app.EnhanceButton.Visible = 'on';
            catch 
                warndlg('Corrupted Image! Please rebrowse.','Warning');
            end
        end

        % Button pushed function: EnhanceButton
        function EnhanceButtonPushed(app, event)
           global ImageFile;
           
            [row,col]=size(ImageFile);
            Output=uint8(zeros(row,col));
            
            kernelValue=app.KernelSizeSpinner.Value;
            kernel=(ones(kernelValue))/(kernelValue*kernelValue);
            [krow,kcol]=size(kernel);
            
            if(size(kernel)>size(ImageFile))
                warndlg('kernel size greater than image size','Warning');
                return;
            end
            %kernel midpoints
            kerMidRow = ceil(krow/2); 
            kerMidCol = ceil(kcol/2);
            
            %padding needed
            zeropad=(kernelValue-1)/2;
            Padding=zeros(row+krow-1,col+kcol-1);
            for i=1:row
                for j=1:col
                        Padding(i+zeropad,j+zeropad)=ImageFile(i,j);
                end
            end
            
            %convolution
            for NewImgRow = kerMidRow : size(Padding,1) - (kerMidRow - 1)
                 for NewImgCol = kerMidCol : size(Padding,2) - (kerMidCol - 1)
                     s = 0;
                     for kerRow = 1 : krow
                         for kerCol = 1 : kcol
                             multVal = kernel(kerRow, kerCol) * Padding((kerRow - kerMidRow + NewImgRow), (kerCol - kerMidCol + NewImgCol)); 
                             s = s + multVal;
                         end
                     end
                     
                     Output((NewImgRow - kerMidRow + 1), (NewImgCol - kerMidCol + 1)) = s;
                 end
            end
            app.BlurredNoisyImageLabel.Text = 'Blurred Image';
            imshow(Output,'Parent' ,app.Blurred_NoisyDisplay);
            
            %enhanced image
            kFactor=app.FactorKEditField.Value;
            enhancedImage=(kFactor*(ImageFile-Output))+ImageFile;
            for i=1:size(enhancedImage,1)
                for j=1:size(enhancedImage,2)

                   if enhancedImage(i,j) > 255
                       enhancedImage(i,j)=255;
                   end
                   if enhancedImage(i,j) < 0
                       enhancedImage(i,j)=0;
                   end 

                end
            end
            app.EnhancedImageLabel.Text = 'Enhanced Image';
            imshow(enhancedImage,'Parent' ,app.EnhancedImageDisplay);

            
        end

        % Button pushed function: RemoveNoiseButton
        function RemoveNoiseButtonPushed(app, event)
            global ImageFile;
            
            [row,col]=size(ImageFile);
            
            % add noise 
            for i=1:row
                for j=1:col
                 NoisyImage(i,j)=ImageFile(i,j);
                end
            end
            
            
            for loop=1:(row*col*0.03)
                i=(ceil(rand*row));
                j=(ceil(rand*col));
                NoisyImage(i,j)=255;
            end
            for loop=1:(row*col*0.03)
                i=(ceil(rand*row));
                j=(ceil(rand*col));
                NoisyImage(i,j)=0;
            end
            app.BlurredNoisyImageLabel.Text = 'Noisy Image';
            imshow(NoisyImage,'Parent' ,app.Blurred_NoisyDisplay);
            
            % remove noise
            NoiseRemoved=uint8(zeros(row,col));
            
            kernelValue=app.KernelSizeSpinner.Value;
            kernel=ones(kernelValue);
            [krow,kcol]=size(kernel);
            
            if(size(kernel)>size(NoisyImage))
                warndlg('kernel size greater than image size','Warning');
                return;
            end
            
            kerMidRow = ceil(krow/2); 
            kerMidCol = ceil(kcol/2);
            
            zeropad=(kernelValue-1)/2;
            Padding=zeros(row+krow-1,col+kcol-1);
            for i=1:row
                for j=1:col
                        Padding(i+zeropad,j+zeropad)=NoisyImage(i,j);
                end
            end
            
            for NewImgRow = kerMidRow : size(Padding,1) - (kerMidRow - 1)
                 for NewImgCol = kerMidCol : size(Padding,2) - (kerMidCol - 1)
                     element=zeros(kernelValue*kernelValue,1);
                     for kerRow = 1 : krow
                         for kerCol = 1 : kcol
                             element(krow*(kerCol-1)+kerRow)=Padding((kerRow - kerMidRow + NewImgRow), (kerCol - kerMidCol + NewImgCol));
                         end
                     end
                     %%%%%%%%%%hangib el median hna
                     sortedValues = zeros(size(element));          
                     for k = 1:numel(element)           
                        [m, ind] = min(element);       
                        element(ind)   = [];            
                        sortedValues(k)= m;
                     end
                     median=sortedValues(ceil(numel(sortedValues)/2));
                     NoiseRemoved((NewImgRow - kerMidRow + 1), (NewImgCol - kerMidCol + 1)) = median;
                 end
            end
            app.EnhancedImageLabel.Text = 'Denoised Image';
            imshow(NoiseRemoved,'Parent' ,app.EnhancedImageDisplay);
        end

        % Button pushed function: FourierBrowseButton
        function FourierBrowseButtonPushed(app, event)
           global fourierWithoutLog;
           global fourierWithLog;
         
            [Filename, Pathname,indx] = uigetfile({'*.jpg';'*.jpeg';'*.bmp';'*.dcm';'*.tif';'*.png'}, 'All Files');
               
            if Filename==0
                return;    
            end
            
            fullname=[Pathname,Filename];
            if (indx==4)
                dcmfile = dicomread(fullname);
                dcmImage = uint8(255 * mat2gray(dcmfile));
                folder='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task4-MariamOkasha/Images';
                imwrite(dcmImage,fullfile(folder,'dcmImage.jpg'));
                fullname='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task4-MariamOkasha/Images/dcmImage.jpg';
            end
       
             try
                ImageFile=imread(fullname);
                [row, col, channels] = size(ImageFile);
                  
                if channels==3
                    ImageFile= rgb2gray(ImageFile);
                end
                imshow(ImageFile,'Parent' ,app.OriginalImageDisplay);
                
                
                % Apply Fourier transform to the image
                fourierWithoutLog=1+fftshift(fft2(ImageFile));
                fourierWithLog=log(fourierWithoutLog);
            
                %get magnitude  
                magnitudeW=abs(fourierWithLog);
                for i=1:size(magnitudeW,1)
                    for j=1:size(magnitudeW,2)
                        
                       if magnitudeW(i,j) > 255
                           magnitudeW(i,j)=255;
                       end
                       if magnitudeW(i,j) < 0
                           magnitudeW(i,j)=0;
                       end 
                       
                    end
                end
                imshow(magnitudeW,[],'Parent' ,app.MagnitudeDisplay);
                app.MagnitudeLabel.Text = 'Magnitude With Log';
                
                %get phase
                phaseW=angle(fourierWithLog);
                for i=1:size(phaseW,1)
                    for j=1:size(phaseW,2)
                        
                       if phaseW(i,j) > 255
                           phaseW(i,j)=255;
                       end
                       if phaseW(i,j) < 0
                           phaseW(i,j)=0;
                       end 
                       
                    end
                end
                imshow(phaseW,[],'Parent' ,app.PhaseDisplay);    
                app.PhaseLabel.Text = 'Phase With Log';        
            catch 
                warndlg('Corrupted Image! Please rebrowse.','Warning');
            end             
        end

        % Value changed function: LogDropDown
        function LogDropDownValueChanged(app, event)
           global fourierWithoutLog;
           global fourierWithLog;
            value = app.LogDropDown.Value;
            
            
            if (strcmp(value,'With Log')==1)
                %get magnitude
                magnitudeW=abs(fourierWithLog);
                for i=1:size(magnitudeW,1)
                    for j=1:size(magnitudeW,2)
                        
                       if magnitudeW(i,j) > 255
                           magnitudeW(i,j)=255;
                       end
                       if magnitudeW(i,j) < 0
                           magnitudeW(i,j)=0;
                       end 
                       
                    end
                end
                imshow(magnitudeW,[],'Parent' ,app.MagnitudeDisplay);
                app.MagnitudeLabel.Text = 'Magnitude With Log';
                %get phase
                phaseW=angle(fourierWithLog);
                for i=1:size(phaseW,1)
                    for j=1:size(phaseW,2)
                        
                       if phaseW(i,j) > 255
                           phaseW(i,j)=255;
                       end
                       if phaseW(i,j) < 0
                           phaseW(i,j)=0;
                       end 
                       
                    end
                end
                imshow(phaseW,[],'Parent' ,app.PhaseDisplay);
                app.PhaseLabel.Text = 'Phase With Log';  
                
            elseif (strcmp(value,'Without Log')==1)
                %get magnitude
                magnitude=abs(fourierWithoutLog);
                for i=1:size(magnitude,1)
                    for j=1:size(magnitude,2)
                        
                       if magnitude(i,j) > 255
                           magnitude(i,j)=255;
                       end
                       if magnitude(i,j) < 0
                           magnitude(i,j)=0;
                       end 
                       
                    end
                end
                imshow(magnitude,[],'Parent' ,app.MagnitudeDisplay);
                app.MagnitudeLabel.Text = 'Magnitude Without Log';
                %get phase
                phase=angle(fourierWithoutLog);
                for i=1:size(phase,1)
                    for j=1:size(phase,2)
                        
                       if phase(i,j) > 255
                           phase(i,j)=255;
                       end
                       if phase(i,j) < 0
                           phase(i,j)=0;
                       end 
                       
                    end
                end
                imshow(phase,[],'Parent' ,app.PhaseDisplay);
                app.PhaseLabel.Text = 'Phase Without Log';
            end
        end

        % Button pushed function: FourierIIBrowseButton
        function FourierIIBrowseButtonPushed(app, event)
          global ImageFileF;
            [Filename, Pathname,indx] = uigetfile({'*.jpg';'*.jpeg';'*.bmp';'*.dcm';'*.tif';'*.png'}, 'All Files');
             
            if Filename==0
                return;    
            end
            
            fullname=[Pathname,Filename];
            if (indx==4)
                dcmfile = dicomread(fullname);
                dcmImage = uint8(255 * mat2gray(dcmfile));
                folder='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task4-MariamOkasha/Images';
                imwrite(dcmImage,fullfile(folder,'dcmImage.jpg'));
                fullname='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task4-MariamOkasha/Images/dcmImage.jpg';
            end
       
            try
                ImageFileF=imread(fullname);
                [row, col, channels] = size(ImageFileF);
                  
                if channels==3
                    ImageFileF= rgb2gray(ImageFileF);
                end
             imshow(ImageFileF,'Parent' ,app.ImageDisplaay);
             
           %padding even images 
            if(rem(col,2)==0)
                ImageFileF(:,end+1)=ImageFileF(:,end);
            end
            if(rem(row,2)==0)
                ImageFileF(end+1,:)=ImageFileF(end,:);
            end 
            [row,col]=size(ImageFileF);
            
            catch 
                warndlg('Corrupted Image! Please rebrowse.','Warning');
            end            
        end

        % Button pushed function: FilterButton
        function FilterButtonPushed(app, event)
            global ImageFileF;
            
            [row,col]=size(ImageFileF);
            kernelValue=app.FourierKernelSizeSpinner.Value;
            kernel=ones(kernelValue)/(kernelValue*kernelValue);
            [krow,kcol]=size(kernel);
            
            if(size(kernel)>size(ImageFileF))
                warndlg('kernel size greater than image size','Warning');
                return;
            end
            Padding=zeros(row,col);
            for i=1:kernelValue
                for j=1:kernelValue
                        Padding(i+round((row-kernelValue)/2),j+round((col-kernelValue)/2))=kernel(i,j);
                end
            end
            
            FourierImage=fftshift((fft2(ImageFileF)));
            FourierFilter=fftshift((fft2(Padding)));
            Mult=FourierImage.*FourierFilter;
            Fourier=uint8(abs(ifftshift(ifft2(Mult))));
            imshow(Fourier,'Parent' ,app.FourierFilterDisplay);
            
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Output=uint8(zeros(row,col));
            kerMidRow = ceil(krow/2); 
            kerMidCol = ceil(kcol/2);
            
            %padding needed
            zeropad=(kernelValue-1)/2;
            PaddingSpatial=zeros(row+krow-1,col+kcol-1);
            for i=1:row
                for j=1:col
                        PaddingSpatial(i+zeropad,j+zeropad)=ImageFileF(i,j);
                end
            end
            %convolution
            for NewImgRow = kerMidRow : size(PaddingSpatial,1) - (kerMidRow - 1)
                 for NewImgCol = kerMidCol : size(PaddingSpatial,2) - (kerMidCol - 1)
                     s = 0;
                     for kerRow = 1 : krow
                         for kerCol = 1 : kcol
                             multVal = kernel(kerRow, kerCol) * PaddingSpatial((kerRow - kerMidRow + NewImgRow), (kerCol - kerMidCol + NewImgCol)); 
                             s = s + multVal;
                         end
                     end
                     
                     Output((NewImgRow - kerMidRow + 1), (NewImgCol - kerMidCol + 1)) = s;
                 end
            end   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           sub=Fourier-Output 
           imshow(sub,'Parent' ,app.DifferenceDisplay);
        end

        % Callback function
        function DifferenceButtonPushed(app, event)
           global ImageFile;
           
            [row,col]=size(ImageFile);
            Output=uint8(zeros(row,col));
            
            kernelValue=app.KernelSizeSpinner.Value;
            kernel=(ones(kernelValue))/(kernelValue*kernelValue);
            [krow,kcol]=size(kernel);
            kerMidRow = ceil(krow/2); 
            kerMidCol = ceil(kcol/2);
            
            %padding needed
            zeropad=(kernelValue-1)/2;
            Padding=zeros(row+krow-1,col+kcol-1);
            for i=1:row
                for j=1:col
                        Padding(i+zeropad,j+zeropad)=ImageFile(i,j);
                end
            end
            
            %convolution
            for NewImgRow = kerMidRow : size(Padding,1) - (kerMidRow - 1)
                 for NewImgCol = kerMidCol : size(Padding,2) - (kerMidCol - 1)
                     s = 0;
                     for kerRow = 1 : krow
                         for kerCol = 1 : kcol
                             multVal = kernel(kerRow, kerCol) * Padding((kerRow - kerMidRow + NewImgRow), (kerCol - kerMidCol + NewImgCol)); 
                             s = s + multVal;
                         end
                     end
                     
                     Output((NewImgRow - kerMidRow + 1), (NewImgCol - kerMidCol + 1)) = s;
                 end
            end
        end

        % Button pushed function: RemovePatternButton
        function RemovePatternButtonPushed(app, event)
            fullname='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task6.2-MariamOkasha/Pattern.jpeg';
            
            ImageFile=imread(fullname);
            [row, col, channels] = size(ImageFile);
            if channels==3
                ImageFile= rgb2gray(ImageFile);
            end
            imshow(ImageFile,'Parent' ,app.ImageWithPattern);
            
            frequencyImage = fftshift(fft2(ImageFile));
            % Take log magnitude so we can see it better in the display.
  


            frequencyImage(482:497,320:335 ) = 0;
            frequencyImage(480:495,375:390) = 0;
            frequencyImage(564:579,322:337 ) = 0;
            frequencyImage(559:574,375:390) = 0; 
            
            
            frequencyImage(447:462,320:335 ) = 0;
            frequencyImage(445:460,375:390) = 0;
            frequencyImage(599:614,322:337 ) = 0;
            frequencyImage(594:609,375:390) = 0; 

            amplitudeImage = log(abs(frequencyImage));
            filteredImage = ifft2(fftshift(frequencyImage));
            amplitudeImage2 = abs(filteredImage);

    
           
              
            imshow(amplitudeImage,[],'Parent' ,app.Mask);
            imshow(amplitudeImage2, [],'Parent' ,app.ImageWithoutPattern);
                      
        end

        % Button pushed function: CreateImageButton
        function CreateImageButtonPushed(app, event)
           global Image;
            Image = 50*ones(256,256, 'uint8');
            Image(35:222,35:222)=150;
            [X,Y] = meshgrid(1:256,1:256) ;
            circle=(X-128.5).^2+(Y-128.5).^2<=60.^2;
            Image(circle)=250;
            app.CreatedImage.Position=[35,85,291,341];
          
            imshow(Image,'Parent' ,app.CreatedImage);
            
            disableDefaultInteractivity(app.CreatedImage)
        end

        % Value changed function: NoiseTypeDropDown
        function NoiseTypeDropDownValueChanged(app, event)
            global OriginalImage;
            global Image;
            global roi;
            global Gaussian_noise;
            global Uniform_noise;
            value = app.NoiseTypeDropDown.Value;
            if (strcmp(value,'No Noise')==1)
                OriginalImage=uint8(Image);
                imshow(OriginalImage,'Parent' ,app.CreatedImage);                 
%                 roi=drawrectangle(app.CreatedImage);
                
            elseif (strcmp(value,'Gaussian Noise')==1)
                
                Gaussian_noise=uint8((sqrt(5^2)*randn(size(Image)))+0+double(Image));
                imshow(Gaussian_noise,'Parent' ,app.CreatedImage); 
%                 roi=drawrectangle(app.CreatedImage);
                
            elseif (strcmp(value,'Uniform Noise')==1)
                a = -10;
                b = 10; 
                Uniform_noise=uint8(a+(b-a)*rand(size(Image))+double(Image));
                imshow(Uniform_noise,'Parent' ,app.CreatedImage); 
%                 roi=drawrectangle(app.CreatedImage);
            end
            roi=drawrectangle(app.CreatedImage);
        end

        % Button pushed function: ConfirmRegionofInterestButton
        function ConfirmRegionofInterestButtonPushed(app, event)
            global roi;
            global OriginalImage;
            global Gaussian_noise;
            global Uniform_noise;
              try
                value = app.NoiseTypeDropDown.Value;
                roi.InteractionsAllowed='none';
                roi_xy=roi.Position;
                x1=round(roi_xy(1));
                x2=round(x1+roi_xy(3));
                y1=round(roi_xy(2));
                y2=round(y1+roi_xy(4)); 
                
                if (strcmp(value,'No Noise')==1)
                    region=OriginalImage(y1:y2,x1:x2);
                    disp(region);
                elseif (strcmp(value,'Gaussian Noise')==1)
                    region=Gaussian_noise(y1:y2,x1:x2); 
                    disp(region);
                elseif (strcmp(value,'Uniform Noise')==1)
                    region=Uniform_noise(y1:y2,x1:x2);
                    
                end
                histo=zeros(1,256);
                [j,k]=size(region);
                for row=1:j
                    for col=1:k
                        
                        pix=region(row,col);
                        histo(pix)=histo(pix)+1;
                        normalized(pix)=histo(pix)/(k*j);
                
                    end
                end 
                 
                s=0;count=0;sumexp=0;
                bar(app.ROIHistogram,normalized);
                for i=1:length(histo)
                    s=s+i*histo(i);
                    sumexp=sumexp+i^2*histo(i);
                    count=count+histo(i);
                end
                mean=s/count;
                Exp=sumexp/count;
                variance=sqrt(Exp-mean^2);
                app.MeanEditField.Value=num2str(mean);
                app.VarianceEditField.Value=num2str(variance);
             catch
                 warndlg('ROI must be inside the Image');
             end
        end

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
            global lam;
            global I2;
            global I3;
            global I4;
            
            Pt = phantom(256);
            imshow(Pt,[],'Parent',app.ScheppLoganPhantom);
            
            theta=0:179;
            an=[0,20,40,60,160];
            
            sino = radon(Pt, theta);
            imshow(rot90(sino),[],'Parent',app.Sinogram);
            
            lam=iradon(radon(Pt,an),an,'linear','none');
            imshow(lam,[],'Parent',app.Laminogram);
            I2 = iradon(sino,theta,'linear','none');
            I3 = iradon(sino, theta, 'linear', 'Ram-Lak');
            I4 = iradon(sino, theta, 'linear', 'Hamming');
            
        end

        % Value changed function: LaminogramDropDown
        function LaminogramDropDownValueChanged(app, event)
            global lam;
            global I2;
            global I3;
            global I4;            
            value = app.LaminogramDropDown.Value;
          
            if (strcmp(value,'Unfiltered Laminogram (0,20,40,60,160)')==1)

                imshow(lam,[],'Parent',app.Laminogram); 
                
            elseif (strcmp(value,'UnFiltered Laminogram (0:179)')==1)
                
                imshow(I2,[],'Parent',app.Laminogram);
                
            elseif (strcmp(value,'RamLak')==1)
                
                imshow(I3,[],'Parent',app.Laminogram);
                
            elseif (strcmp(value,'Hamming')==1)
                
                imshow(I4,[],'Parent',app.Laminogram);
            end
        end

        % Button pushed function: ManualButton
        function ManualButtonPushed(app, event)
            Pt = phantom(256);
            [x,y] = size(Pt);
            %Set a matrix R to hold the projection data
            R = zeros((180),x);
            %Set vector for the angles to take projections at
            angles = 0:179;
            %Matrix to hold all the rotated images separately to prevent blurring
            rotatedImage = zeros(x,y,length(angles));
        
            %Loop to rotate image and add up values for the projections
            for i = 1:(length(angles))
                %rotate the image starting with theta = 0 degrees
                rotatedImage(:,:,i) = imrotate(Pt,angles(i),'nearest','crop');
                %Sum the columns of img to get projection data
                %Each row of R contains a projection at a certain theta
                R(i,:) = sum(rotatedImage(:,:,i),1);
            end

            %Convert the matrix to a gray scale image in the range 0 to 255
            R = mat2gray(R);

            imshow(flip(R,2),[],'Parent',app.Sinogram);
            

        end

        % Button pushed function: BuiltinButton
        function BuiltinButtonPushed(app, event)
            Pt = phantom(256);
            sino = radon(Pt);
            imshow(rot90(sino),[],'Parent',app.Sinogram);            
        end

        % Button pushed function: DisplayButton
        function DisplayButtonPushed(app, event)
            global StructElement;
            global ImageFile;
            fullname='E:/CUFE/Fall_2022/Medical Image Processing & Computer Vision/Tasks/Task9-MariamOkasha/binary_image.png';
            ImageFile=imread(fullname);
            ImageFile= imbinarize(rgb2gray(ImageFile));
            imshow(ImageFile,'Parent' ,app.FingerprintImage);
            StructElement = zeros(5,5);
            StructElement(2:4, 1:5) = 1;
            StructElement(1:5, 2:4) = 1;
            imshow(StructElement,'Parent' ,app.MorphologicalImage);  
    
        end

        % Button pushed function: ErosionButton
        function ErosionButtonPushed(app, event)
            global StructElement;
            global ImageFile;
            [row,col]=size(ImageFile);  
            imgErode= zeros(row,col);
            
            for i=3:row-2
                for j=3:col-2
                    s=0;
                    temp= ImageFile(i-2:i+2, j-2:j+2);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==21)
                        imgErode(i,j)= 1;
                    end
                    
                end
            end
            imshow((imgErode),'Parent' ,app.MorphologicalImage); 
        end

        % Button pushed function: DilationButton
        function DilationButtonPushed(app, event)
            global StructElement;
            global ImageFile;
            [row,col]=size(ImageFile);  
            imgDilate= zeros(row,col);
            Padding=zeros(row+4,col+4);
            [newrow,newcol]=size(Padding);
            for i=1:row
                for j=1:col
                        Padding(i+2,j+2)=ImageFile(i,j);
                end
            end
            
            for i=3:newrow-2
                for j=3:newcol-2
                    temp= Padding(i-2:i+2, j-2:j+2);
                    product= temp.*StructElement;
                    imgDilate(i-2,j-2)= max(max(product));
                end
            end
            
            imshow((imgDilate),'Parent' ,app.MorphologicalImage);     
        end

        % Button pushed function: OpeningButton
        function OpeningButtonPushed(app, event)
            global StructElement;
            global ImageFile;
            [row,col]=size(ImageFile);  
            imgErode= zeros(row,col);
            imgOpening= zeros(row,col);
            Padding=zeros(row+4,col+4);
            [newrow,newcol]=size(Padding);
            
            for i=3:row-2
                for j=3:col-2
                    s=0;
                    temp= ImageFile(i-2:i+2, j-2:j+2);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==21)
                        imgErode(i,j)= 1;
                    end
                    
                end
            end
            
            for i=1:row
                for j=1:col
                        Padding(i+2,j+2)=imgErode(i,j);
                end
            end
            
            for i=3:newrow-2
                for j=3:newcol-2
                    temp= Padding(i-2:i+2, j-2:j+2);
                    product= temp.*StructElement;
                    imgOpening(i-2,j-2)= max(max(product));
                end
            end
            
            imshow((imgOpening),'Parent' ,app.MorphologicalImage);             
            
        end

        % Button pushed function: ClosingButton
        function ClosingButtonPushed(app, event)
            global StructElement;
            global ImageFile;
            [row,col]=size(ImageFile);  
            imgDilate= zeros(row,col);
            imgClosing= zeros(row,col);
            Padding=zeros(row+4,col+4);
            [newrow,newcol]=size(Padding);  
            for i=1:row
                for j=1:col
                        Padding(i+2,j+2)=ImageFile(i,j);
                end
            end
            
            for i=3:newrow-2
                for j=3:newcol-2
                    temp= Padding(i-2:i+2, j-2:j+2);
                    product= temp.*StructElement;
                    imgDilate(i-2,j-2)= max(max(product));
                end
            end
            
            for i=3:row-2
                for j=3:col-2
                    s=0;
                    temp= imgDilate(i-2:i+2, j-2:j+2);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==21)
                        imgClosing(i,j)= 1;
                    end
                    
                end
            end  
            imshow((imgClosing),'Parent' ,app.MorphologicalImage);             
            
        end

        % Button pushed function: Button
        function ButtonPushed(app, event)
            global ImageFile;
            StructElement = zeros(3,3);
            StructElement(2, 1:3) = 1;
            StructElement(1:3, 2) = 1; 
            [row,col]=size(ImageFile);  
            imgErode= zeros(row,col);
            imgOpening= zeros(row,col);
            imgDilate= zeros(row,col);
            imgClosing= zeros(row,col);
            Padding1=zeros(row+2,col+2);
            Padding2=zeros(row+2,col+2);
            [newrow,newcol]=size(Padding1);
            
            for i=2:row-1
                for j=2:col-1
                    s=0;
                    temp= ImageFile(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==5)
                        imgErode(i,j)= 1;
                    end
                    
                end
            end
            
            for i=1:row
                for j=1:col
                        Padding1(i+1,j+1)=imgErode(i,j);
                end
            end
            
            for i=2:newrow-1
                for j=2:newcol-1
                    temp= Padding1(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    imgOpening(i-1,j-1)= max(max(product));
                end
            end  
            imshow((imgOpening),'Parent' ,app.MorphologicalImage);
            
            for i=1:row
                for j=1:col
                        Padding2(i+1,j+1)=imgOpening(i,j);
                end
            end
            
            for i=2:newrow-1
                for j=2:newcol-1
                    temp= Padding2(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    imgDilate(i-1,j-1)= max(max(product));
                end
            end 
            
            for i=2:row-1
                for j=2:col-1
                    s=0;
                    temp= imgDilate(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==5)
                        imgClosing(i,j)= 1;
                    end
                    
                end
            end 
            imshow((imgClosing),'Parent' ,app.MorphologicalImage);
        end

        % Button pushed function: Button_2
        function Button_2Pushed(app, event)
            global ImageFile;
            StructElement = ones(1,3);
    
            [row,col]=size(ImageFile);  
            imgErode= zeros(row,col);
            imgOpening= zeros(row,col);
            imgDilate= zeros(row,col);
            imgClosing= zeros(row,col);
            Padding1=zeros(row,col+2);
            Padding2=zeros(row,col+2);
            [newrow,newcol]=size(Padding1);
            
            for i=1:row
                for j=2:col-1
                    s=0;
                    temp= ImageFile(i:i, j-1:j+1);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==3)
                        imgErode(i,j)= 1;
                    end
                    
                end
            end
            
            for i=1:row
                for j=1:col
                        Padding1(i,j+1)=imgErode(i,j);
                end
            end
            
            for i=1:newrow
                for j=2:newcol-1
                    temp= Padding1(i:i, j-1:j+1);
                    product= temp.*StructElement;
                    imgOpening(i,j-1)= max(max(product));
                end
            end  
            imshow((imgOpening),'Parent' ,app.MorphologicalImage);
            
            for i=1:row
                for j=1:col
                        Padding2(i,j+1)=imgOpening(i,j);
                end
            end
            
            for i=1:newrow
                for j=2:newcol-1
                    temp= Padding2(i:i, j-1:j+1);
                    product= temp.*StructElement;
                    imgDilate(i,j-1)= max(max(product));
                end
            end 
            
            for i=1:row
                for j=2:col-1
                    s=0;
                    temp= imgDilate(i:i, j-1:j+1);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==3)
                        imgClosing(i,j)= 1;
                    end
                    
                end
            end 
            imshow((imgClosing),'Parent' ,app.MorphologicalImage);       
        end

        % Button pushed function: Button_3
        function Button_3Pushed(app, event)
            global ImageFile;
            global imgClosinggg;
            StructElement = ones(3,3);
            StructElement(1, 1) = 0;
            StructElement(3, 3) = 0; 
            [row,col]=size(ImageFile);  
            imgErode= zeros(row,col);
            imgOpening= zeros(row,col);
            imgDilate= zeros(row,col);
            imgClosinggg= zeros(row,col);
            Padding1=zeros(row+2,col+2);
            Padding2=zeros(row+2,col+2);
            [newrow,newcol]=size(Padding1);
            
            for i=2:row-1
                for j=2:col-1
                    s=0;
                    temp= ImageFile(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==7)
                        imgErode(i,j)= 1;
                    end
                    
                end
            end
            
            for i=1:row
                for j=1:col
                        Padding1(i+1,j+1)=imgErode(i,j);
                end
            end
            
            for i=2:newrow-1
                for j=2:newcol-1
                    temp= Padding1(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    imgOpening(i-1,j-1)= max(max(product));
                end
            end  
            imshow((imgOpening),'Parent' ,app.MorphologicalImage);
            
            for i=1:row
                for j=1:col
                        Padding2(i+1,j+1)=imgOpening(i,j);
                end
            end
            
            for i=2:newrow-1
                for j=2:newcol-1
                    temp= Padding2(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    imgDilate(i-1,j-1)= max(max(product));
                end
            end 
            
            for i=2:row-1
                for j=2:col-1
                    s=0;
                    temp= imgDilate(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==7)
                        imgClosinggg(i,j)= 1;
                    end
                    
                end
            end 
            imshow((imgClosinggg),'Parent' ,app.MorphologicalImage);            
        end

        % Button pushed function: Button_4
        function Button_4Pushed(app, event)
            global ImageFile;
            global imgClosingg;
            StructElement = ones(3,3);
            StructElement(3, 1) = 0;
            StructElement(1, 3) = 0; 
            [row,col]=size(ImageFile);  
            imgErode= zeros(row,col);
            imgOpening= zeros(row,col);
            imgDilate= zeros(row,col);
            imgClosingg= zeros(row,col);
            Padding1=zeros(row+2,col+2);
            Padding2=zeros(row+2,col+2);
            [newrow,newcol]=size(Padding1);
            
            for i=2:row-1
                for j=2:col-1
                    s=0;
                    temp= ImageFile(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==7)
                        imgErode(i,j)= 1;
                    end
                    
                end
            end
            
            for i=1:row
                for j=1:col
                        Padding1(i+1,j+1)=imgErode(i,j);
                end
            end
            
            for i=2:newrow-1
                for j=2:newcol-1
                    temp= Padding1(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    imgOpening(i-1,j-1)= max(max(product));
                end
            end  
            imshow((imgOpening),'Parent' ,app.MorphologicalImage);
            
            for i=1:row
                for j=1:col
                        Padding2(i+1,j+1)=imgOpening(i,j);
                end
            end
            
            for i=2:newrow-1
                for j=2:newcol-1
                    temp= Padding2(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    imgDilate(i-1,j-1)= max(max(product));
                end
            end 
            
            for i=2:row-1
                for j=2:col-1
                    s=0;
                    temp= imgDilate(i-1:i+1, j-1:j+1);
                    product= temp.*StructElement;
                    s=sum(product,"all");
                    if (s==7)
                        imgClosingg(i,j)= 1;
                    end
                    
                end
            end 
            imshow((imgClosingg),'Parent' ,app.MorphologicalImage);
        end

        % Callback function
        function Button_5Pushed(app, event)
            global imgClosingg; 
            global imgClosinggg;
            x=max(imgClosingg,imgClosinggg);
            imshow((x),'Parent' ,app.MorphologicalImage);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ImageProcessing and hide until all components are created
            app.ImageProcessing = uifigure('Visible', 'off');
            app.ImageProcessing.Color = [0.9412 0.9412 0.9412];
            app.ImageProcessing.Position = [100 100 1016 560];
            app.ImageProcessing.Name = 'Image Processing';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.ImageProcessing);
            app.TabGroup.Position = [1 3 1016 559];

            % Create ImageViewerTab
            app.ImageViewerTab = uitab(app.TabGroup);
            app.ImageViewerTab.Title = 'Image Viewer';
            app.ImageViewerTab.BackgroundColor = [0.902 0.949 0.8902];

            % Create BrowseButton
            app.BrowseButton = uibutton(app.ImageViewerTab, 'push');
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);
            app.BrowseButton.FontName = 'Yu Gothic Medium';
            app.BrowseButton.FontSize = 16;
            app.BrowseButton.Position = [44 476 89 34];
            app.BrowseButton.Text = 'Browse';

            % Create PathEditField
            app.PathEditField = uieditfield(app.ImageViewerTab, 'text');
            app.PathEditField.FontName = 'Helvetica-Narrow';
            app.PathEditField.Position = [44 416 616 26];

            % Create WidthLabel
            app.WidthLabel = uilabel(app.ImageViewerTab);
            app.WidthLabel.FontName = 'Yu Gothic';
            app.WidthLabel.FontSize = 16;
            app.WidthLabel.Position = [463 285 89 32];
            app.WidthLabel.Text = 'Width';

            % Create HeightLabel
            app.HeightLabel = uilabel(app.ImageViewerTab);
            app.HeightLabel.FontName = 'Yu Gothic';
            app.HeightLabel.FontSize = 16;
            app.HeightLabel.Position = [463 244 89 28];
            app.HeightLabel.Text = 'Height';

            % Create SizeinbitsLabel
            app.SizeinbitsLabel = uilabel(app.ImageViewerTab);
            app.SizeinbitsLabel.FontName = 'Yu Gothic';
            app.SizeinbitsLabel.FontSize = 16;
            app.SizeinbitsLabel.Position = [465 145 85 26];
            app.SizeinbitsLabel.Text = 'Size in bits';

            % Create BitDepthLabel
            app.BitDepthLabel = uilabel(app.ImageViewerTab);
            app.BitDepthLabel.FontName = 'Yu Gothic';
            app.BitDepthLabel.FontSize = 16;
            app.BitDepthLabel.Position = [465 194 89 28];
            app.BitDepthLabel.Text = 'Bit Depth';

            % Create ImageColorLabel
            app.ImageColorLabel = uilabel(app.ImageViewerTab);
            app.ImageColorLabel.FontName = 'Yu Gothic';
            app.ImageColorLabel.FontSize = 16;
            app.ImageColorLabel.Position = [465 98 94 26];
            app.ImageColorLabel.Text = 'Image Color';

            % Create WidthEditField
            app.WidthEditField = uieditfield(app.ImageViewerTab, 'text');
            app.WidthEditField.FontName = 'Yu Gothic';
            app.WidthEditField.FontSize = 14;
            app.WidthEditField.Position = [576 293 80 27];

            % Create HeightEditField
            app.HeightEditField = uieditfield(app.ImageViewerTab, 'text');
            app.HeightEditField.FontName = 'Yu Gothic';
            app.HeightEditField.FontSize = 14;
            app.HeightEditField.Position = [576 244 80 27];

            % Create SizeInBitsEditField
            app.SizeInBitsEditField = uieditfield(app.ImageViewerTab, 'text');
            app.SizeInBitsEditField.FontName = 'Yu Gothic';
            app.SizeInBitsEditField.FontSize = 14;
            app.SizeInBitsEditField.Position = [576 144 80 27];

            % Create BitDepthEditField
            app.BitDepthEditField = uieditfield(app.ImageViewerTab, 'text');
            app.BitDepthEditField.FontName = 'Yu Gothic';
            app.BitDepthEditField.FontSize = 14;
            app.BitDepthEditField.Position = [576 194 80 27];

            % Create ImageColorEditField
            app.ImageColorEditField = uieditfield(app.ImageViewerTab, 'text');
            app.ImageColorEditField.FontName = 'Yu Gothic';
            app.ImageColorEditField.FontSize = 14;
            app.ImageColorEditField.Position = [576 97 80 27];

            % Create ModalityLabel
            app.ModalityLabel = uilabel(app.ImageViewerTab);
            app.ModalityLabel.FontName = 'Yu Gothic';
            app.ModalityLabel.FontSize = 16;
            app.ModalityLabel.Position = [688 289 89 28];
            app.ModalityLabel.Text = 'Modality';

            % Create PatientNameLabel
            app.PatientNameLabel = uilabel(app.ImageViewerTab);
            app.PatientNameLabel.FontName = 'Yu Gothic';
            app.PatientNameLabel.FontSize = 16;
            app.PatientNameLabel.Position = [688 245 102 26];
            app.PatientNameLabel.Text = 'PatientName';

            % Create PatientAgeLabel
            app.PatientAgeLabel = uilabel(app.ImageViewerTab);
            app.PatientAgeLabel.FontName = 'Yu Gothic';
            app.PatientAgeLabel.FontSize = 16;
            app.PatientAgeLabel.Position = [688 197 86 26];
            app.PatientAgeLabel.Text = 'PatientAge';

            % Create BodyPartExaminedLabel
            app.BodyPartExaminedLabel = uilabel(app.ImageViewerTab);
            app.BodyPartExaminedLabel.FontName = 'Yu Gothic';
            app.BodyPartExaminedLabel.FontSize = 16;
            app.BodyPartExaminedLabel.Position = [688 147 145 26];
            app.BodyPartExaminedLabel.Text = 'BodyPartExamined';

            % Create ModalityEditField
            app.ModalityEditField = uieditfield(app.ImageViewerTab, 'text');
            app.ModalityEditField.FontName = 'Yu Gothic';
            app.ModalityEditField.FontSize = 14;
            app.ModalityEditField.Position = [803 289 171 27];

            % Create PatientNameEditField
            app.PatientNameEditField = uieditfield(app.ImageViewerTab, 'text');
            app.PatientNameEditField.FontName = 'Yu Gothic';
            app.PatientNameEditField.FontSize = 14;
            app.PatientNameEditField.Position = [803 244 171 27];

            % Create PatientAgeEditField
            app.PatientAgeEditField = uieditfield(app.ImageViewerTab, 'text');
            app.PatientAgeEditField.FontName = 'Yu Gothic';
            app.PatientAgeEditField.FontSize = 14;
            app.PatientAgeEditField.Position = [803 194 171 27];

            % Create BodyPartExaminedEditField
            app.BodyPartExaminedEditField = uieditfield(app.ImageViewerTab, 'text');
            app.BodyPartExaminedEditField.FontName = 'Yu Gothic';
            app.BodyPartExaminedEditField.FontSize = 14;
            app.BodyPartExaminedEditField.Position = [844 146 130 27];

            % Create ImageDisplay
            app.ImageDisplay = uiaxes(app.ImageViewerTab);
            zlabel(app.ImageDisplay, 'Z')
            app.ImageDisplay.PlotBoxAspectRatio = [1.5948275862069 1 1];
            app.ImageDisplay.XColor = [1 1 1];
            app.ImageDisplay.XTick = [];
            app.ImageDisplay.YColor = [1 1 1];
            app.ImageDisplay.YTick = [];
            app.ImageDisplay.Position = [39 47 395 339];

            % Create InterpolationTab
            app.InterpolationTab = uitab(app.TabGroup);
            app.InterpolationTab.Title = 'Interpolation';
            app.InterpolationTab.BackgroundColor = [0.8588 0.9294 0.9608];

            % Create InterpolateButton
            app.InterpolateButton = uibutton(app.InterpolationTab, 'push');
            app.InterpolateButton.ButtonPushedFcn = createCallbackFcn(app, @InterpolateButtonPushed, true);
            app.InterpolateButton.FontName = 'Yu Gothic Medium';
            app.InterpolateButton.FontSize = 16;
            app.InterpolateButton.Position = [44 476 100 34];
            app.InterpolateButton.Text = 'Interpolate';

            % Create ZoomingFactorLabel
            app.ZoomingFactorLabel = uilabel(app.InterpolationTab);
            app.ZoomingFactorLabel.FontName = 'Yu Gothic';
            app.ZoomingFactorLabel.FontSize = 16;
            app.ZoomingFactorLabel.Visible = 'off';
            app.ZoomingFactorLabel.Position = [232 477 120 28];
            app.ZoomingFactorLabel.Text = 'Zooming Factor';

            % Create ZoomingFactorEditField
            app.ZoomingFactorEditField = uieditfield(app.InterpolationTab, 'numeric');
            app.ZoomingFactorEditField.Limits = [0 16];
            app.ZoomingFactorEditField.FontName = 'Yu Gothic';
            app.ZoomingFactorEditField.FontSize = 14;
            app.ZoomingFactorEditField.Visible = 'off';
            app.ZoomingFactorEditField.Position = [361 481 47 21];

            % Create InterpolationTypeDropDown
            app.InterpolationTypeDropDown = uidropdown(app.InterpolationTab);
            app.InterpolationTypeDropDown.Items = {'Select Interpolation Method', 'Nearest-Neighbour Interpolation', 'Linear Interpolation'};
            app.InterpolationTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @InterpolationTypeDropDownValueChanged, true);
            app.InterpolationTypeDropDown.Visible = 'off';
            app.InterpolationTypeDropDown.FontName = 'Yu Gothic';
            app.InterpolationTypeDropDown.FontSize = 16;
            app.InterpolationTypeDropDown.BackgroundColor = [0.8196 0.8392 0.851];
            app.InterpolationTypeDropDown.Position = [551 474 439 36];
            app.InterpolationTypeDropDown.Value = 'Select Interpolation Method';

            % Create InterpolationPanel
            app.InterpolationPanel = uipanel(app.InterpolationTab);
            app.InterpolationPanel.BorderType = 'none';
            app.InterpolationPanel.BackgroundColor = [0.8588 0.9294 0.9608];
            app.InterpolationPanel.Scrollable = 'on';
            app.InterpolationPanel.Position = [515 195 459 227];

            % Create InterpolatedImage
            app.InterpolatedImage = uiaxes(app.InterpolationPanel);
            zlabel(app.InterpolatedImage, 'Z')
            app.InterpolatedImage.PlotBoxAspectRatio = [2.21694915254237 1 1];
            app.InterpolatedImage.XColor = [1 1 1];
            app.InterpolatedImage.XTick = [];
            app.InterpolatedImage.YColor = [1 1 1];
            app.InterpolatedImage.YTick = [];
            app.InterpolatedImage.Clipping = 'off';
            app.InterpolatedImage.Position = [23 208 10 10];

            % Create OriginalPanel
            app.OriginalPanel = uipanel(app.InterpolationTab);
            app.OriginalPanel.BorderType = 'none';
            app.OriginalPanel.BackgroundColor = [0.8588 0.9294 0.9608];
            app.OriginalPanel.Scrollable = 'on';
            app.OriginalPanel.Position = [44 195 453 228];

            % Create OriginalImage
            app.OriginalImage = uiaxes(app.OriginalPanel);
            zlabel(app.OriginalImage, 'Z')
            app.OriginalImage.XColor = 'none';
            app.OriginalImage.XTick = [];
            app.OriginalImage.YColor = 'none';
            app.OriginalImage.YTick = [];
            app.OriginalImage.ZColor = [1 1 1];
            app.OriginalImage.Clipping = 'off';
            app.OriginalImage.Position = [14 212 10 10];

            % Create OriginalImageinGrayscaleLabel
            app.OriginalImageinGrayscaleLabel = uilabel(app.InterpolationTab);
            app.OriginalImageinGrayscaleLabel.HorizontalAlignment = 'center';
            app.OriginalImageinGrayscaleLabel.FontName = 'Yu Gothic Medium';
            app.OriginalImageinGrayscaleLabel.FontSize = 20;
            app.OriginalImageinGrayscaleLabel.Position = [58 421 428 35];
            app.OriginalImageinGrayscaleLabel.Text = 'Original Image in Grayscale';

            % Create InterpolatedImageLabel
            app.InterpolatedImageLabel = uilabel(app.InterpolationTab);
            app.InterpolatedImageLabel.HorizontalAlignment = 'center';
            app.InterpolatedImageLabel.FontName = 'Yu Gothic Medium';
            app.InterpolatedImageLabel.FontSize = 20;
            app.InterpolatedImageLabel.Position = [533 421 428 35];
            app.InterpolatedImageLabel.Text = 'Interpolated Image';

            % Create OriginalWidthLabel
            app.OriginalWidthLabel = uilabel(app.InterpolationTab);
            app.OriginalWidthLabel.FontName = 'Yu Gothic';
            app.OriginalWidthLabel.FontSize = 16;
            app.OriginalWidthLabel.Position = [94 125 120 28];
            app.OriginalWidthLabel.Text = 'Original Width';

            % Create OriginalHeightLabel
            app.OriginalHeightLabel = uilabel(app.InterpolationTab);
            app.OriginalHeightLabel.FontName = 'Yu Gothic';
            app.OriginalHeightLabel.FontSize = 16;
            app.OriginalHeightLabel.Position = [94 77 120 28];
            app.OriginalHeightLabel.Text = 'Original Height';

            % Create NewWidthLabel
            app.NewWidthLabel = uilabel(app.InterpolationTab);
            app.NewWidthLabel.FontName = 'Yu Gothic';
            app.NewWidthLabel.FontSize = 16;
            app.NewWidthLabel.Position = [591 125 120 28];
            app.NewWidthLabel.Text = 'New Width';

            % Create NewHeightLabel
            app.NewHeightLabel = uilabel(app.InterpolationTab);
            app.NewHeightLabel.FontName = 'Yu Gothic';
            app.NewHeightLabel.FontSize = 16;
            app.NewHeightLabel.Position = [592 77 120 28];
            app.NewHeightLabel.Text = 'New Height';

            % Create OriginalWidthEditField
            app.OriginalWidthEditField = uieditfield(app.InterpolationTab, 'text');
            app.OriginalWidthEditField.FontName = 'Yu Gothic';
            app.OriginalWidthEditField.FontSize = 14;
            app.OriginalWidthEditField.Position = [275 124 140 27];

            % Create OriginalHeightEditField
            app.OriginalHeightEditField = uieditfield(app.InterpolationTab, 'text');
            app.OriginalHeightEditField.FontName = 'Yu Gothic';
            app.OriginalHeightEditField.FontSize = 14;
            app.OriginalHeightEditField.Position = [275 79 140 27];

            % Create NewWidthEditField
            app.NewWidthEditField = uieditfield(app.InterpolationTab, 'text');
            app.NewWidthEditField.FontName = 'Yu Gothic';
            app.NewWidthEditField.FontSize = 14;
            app.NewWidthEditField.Position = [737 126 153 27];

            % Create NewHeightEditField
            app.NewHeightEditField = uieditfield(app.InterpolationTab, 'text');
            app.NewHeightEditField.FontName = 'Yu Gothic';
            app.NewHeightEditField.FontSize = 14;
            app.NewHeightEditField.Position = [737 78 153 27];

            % Create info
            app.info = uilabel(app.InterpolationTab);
            app.info.FontName = 'Yu Gothic UI';
            app.info.FontSize = 8;
            app.info.FontWeight = 'bold';
            app.info.Position = [94 15 762 22];
            app.info.Text = '*Interpolate Button transfers Photo from Image Viewer to Interpolation tab, once clicked while there is a valid image in image viewer tab, user will be able to choose zooming factor and Interpolation type.';

            % Create AffineTransformationTab
            app.AffineTransformationTab = uitab(app.TabGroup);
            app.AffineTransformationTab.Title = 'Affine Transformation';
            app.AffineTransformationTab.BackgroundColor = [0.949 0.9529 1];

            % Create ConstructImageButton
            app.ConstructImageButton = uibutton(app.AffineTransformationTab, 'push');
            app.ConstructImageButton.ButtonPushedFcn = createCallbackFcn(app, @ConstructImageButtonPushed, true);
            app.ConstructImageButton.FontName = 'Yu Gothic Medium';
            app.ConstructImageButton.FontSize = 16;
            app.ConstructImageButton.Position = [43 476 137 34];
            app.ConstructImageButton.Text = 'Construct Image';

            % Create AngleSlider
            app.AngleSlider = uislider(app.AffineTransformationTab);
            app.AngleSlider.Limits = [-180 180];
            app.AngleSlider.MajorTicks = [-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180];
            app.AngleSlider.MajorTickLabels = {'-180', '-150', '-120', '-90', '-60', '-30', '0', '30', '60', '90', '120', '150', '180'};
            app.AngleSlider.ValueChangedFcn = createCallbackFcn(app, @AngleSliderValueChanged, true);
            app.AngleSlider.MinorTicks = [-180 -175 -170 -165 -160 -155 -150 -145 -140 -135 -130 -125 -120 -115 -110 -105 -100 -95 -90 -85 -80 -75 -70 -65 -60 -55 -50 -45 -40 -35 -30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130 135 140 145 150 155 160 165 170 175 180];
            app.AngleSlider.FontName = 'Yu Gothic';
            app.AngleSlider.Position = [567 262 383 7];

            % Create AngleSliderLabel
            app.AngleSliderLabel = uilabel(app.AffineTransformationTab);
            app.AngleSliderLabel.FontName = 'Yu Gothic';
            app.AngleSliderLabel.FontSize = 16;
            app.AngleSliderLabel.Position = [403 241 120 28];
            app.AngleSliderLabel.Text = 'Angle Slider';

            % Create InterpolationTypeOfRotatedImageDropDown
            app.InterpolationTypeOfRotatedImageDropDown = uidropdown(app.AffineTransformationTab);
            app.InterpolationTypeOfRotatedImageDropDown.Items = {'Select Interpolation Method', 'Nearest-Neighbour Interpolation', 'Linear Interpolation'};
            app.InterpolationTypeOfRotatedImageDropDown.ValueChangedFcn = createCallbackFcn(app, @InterpolationTypeOfRotatedImageDropDownValueChanged, true);
            app.InterpolationTypeOfRotatedImageDropDown.FontName = 'Yu Gothic';
            app.InterpolationTypeOfRotatedImageDropDown.FontSize = 14;
            app.InterpolationTypeOfRotatedImageDropDown.BackgroundColor = [0.8196 0.8392 0.851];
            app.InterpolationTypeOfRotatedImageDropDown.Position = [567 149 383 32];
            app.InterpolationTypeOfRotatedImageDropDown.Value = 'Select Interpolation Method';

            % Create InterpolationTypeLabel
            app.InterpolationTypeLabel = uilabel(app.AffineTransformationTab);
            app.InterpolationTypeLabel.FontName = 'Yu Gothic';
            app.InterpolationTypeLabel.FontSize = 16;
            app.InterpolationTypeLabel.Position = [403 152 138 28];
            app.InterpolationTypeLabel.Text = 'Interpolation Type';

            % Create ShearingAngleLabel
            app.ShearingAngleLabel = uilabel(app.AffineTransformationTab);
            app.ShearingAngleLabel.FontName = 'Yu Gothic';
            app.ShearingAngleLabel.FontSize = 16;
            app.ShearingAngleLabel.Position = [403 71 116 28];
            app.ShearingAngleLabel.Text = 'Shearing Angle';

            % Create Shear45Button
            app.Shear45Button = uibutton(app.AffineTransformationTab, 'push');
            app.Shear45Button.ButtonPushedFcn = createCallbackFcn(app, @Shear45ButtonPushed, true);
            app.Shear45Button.BackgroundColor = [0.8196 0.8392 0.851];
            app.Shear45Button.FontName = 'Yu Gothic';
            app.Shear45Button.FontSize = 14;
            app.Shear45Button.Position = [567 66 167 30];
            app.Shear45Button.Text = 'Shear 45';

            % Create ShearNegative45Button
            app.ShearNegative45Button = uibutton(app.AffineTransformationTab, 'push');
            app.ShearNegative45Button.ButtonPushedFcn = createCallbackFcn(app, @ShearNegative45ButtonPushed, true);
            app.ShearNegative45Button.BackgroundColor = [0.8196 0.8392 0.851];
            app.ShearNegative45Button.FontName = 'Yu Gothic';
            app.ShearNegative45Button.FontSize = 14;
            app.ShearNegative45Button.Position = [783 66 167 30];
            app.ShearNegative45Button.Text = 'Shear -45';

            % Create AngleLabel
            app.AngleLabel = uilabel(app.AffineTransformationTab);
            app.AngleLabel.FontName = 'Yu Gothic';
            app.AngleLabel.FontSize = 16;
            app.AngleLabel.Position = [403 375 120 28];
            app.AngleLabel.Text = 'Angle';

            % Create DirectionLabel
            app.DirectionLabel = uilabel(app.AffineTransformationTab);
            app.DirectionLabel.FontName = 'Yu Gothic';
            app.DirectionLabel.FontSize = 16;
            app.DirectionLabel.Position = [401 322 120 28];
            app.DirectionLabel.Text = 'Direction';

            % Create AngleEditField
            app.AngleEditField = uieditfield(app.AffineTransformationTab, 'text');
            app.AngleEditField.FontName = 'Yu Gothic';
            app.AngleEditField.FontSize = 14;
            app.AngleEditField.Position = [503 373 123 27];

            % Create DirectionEditField
            app.DirectionEditField = uieditfield(app.AffineTransformationTab, 'text');
            app.DirectionEditField.FontName = 'Yu Gothic';
            app.DirectionEditField.FontSize = 14;
            app.DirectionEditField.Position = [503 319 123 27];

            % Create RotateWidthLabel
            app.RotateWidthLabel = uilabel(app.AffineTransformationTab);
            app.RotateWidthLabel.FontName = 'Yu Gothic';
            app.RotateWidthLabel.FontSize = 16;
            app.RotateWidthLabel.Position = [719 372 120 28];
            app.RotateWidthLabel.Text = 'Width';

            % Create RotateHeightLabel
            app.RotateHeightLabel = uilabel(app.AffineTransformationTab);
            app.RotateHeightLabel.FontName = 'Yu Gothic';
            app.RotateHeightLabel.FontSize = 16;
            app.RotateHeightLabel.Position = [717 319 120 28];
            app.RotateHeightLabel.Text = 'Height';

            % Create RotateWidthEditField
            app.RotateWidthEditField = uieditfield(app.AffineTransformationTab, 'text');
            app.RotateWidthEditField.FontName = 'Yu Gothic';
            app.RotateWidthEditField.FontSize = 14;
            app.RotateWidthEditField.Position = [819 370 123 27];

            % Create RotateHeightEditField
            app.RotateHeightEditField = uieditfield(app.AffineTransformationTab, 'text');
            app.RotateHeightEditField.FontName = 'Yu Gothic';
            app.RotateHeightEditField.FontSize = 14;
            app.RotateHeightEditField.Position = [819 316 123 27];

            % Create TImage
            app.TImage = uiaxes(app.AffineTransformationTab);
            zlabel(app.TImage, 'Z')
            app.TImage.PlotBoxAspectRatio = [1.06441717791411 1 1];
            app.TImage.YDir = 'reverse';
            app.TImage.XAxisLocation = 'top';
            app.TImage.XColor = [1 1 1];
            app.TImage.XTick = [];
            app.TImage.YColor = [1 1 1];
            app.TImage.YTick = [];
            app.TImage.Position = [35 57 317 385];

            % Create HistogramTab
            app.HistogramTab = uitab(app.TabGroup);
            app.HistogramTab.Title = 'Histogram';
            app.HistogramTab.BackgroundColor = [1 0.9961 0.8314];

            % Create HistogramBrowseButton
            app.HistogramBrowseButton = uibutton(app.HistogramTab, 'push');
            app.HistogramBrowseButton.ButtonPushedFcn = createCallbackFcn(app, @HistogramBrowseButtonPushed, true);
            app.HistogramBrowseButton.FontName = 'Yu Gothic Medium';
            app.HistogramBrowseButton.FontSize = 16;
            app.HistogramBrowseButton.Position = [44 476 89 34];
            app.HistogramBrowseButton.Text = 'Browse';

            % Create EqualizedImageLabel
            app.EqualizedImageLabel = uilabel(app.HistogramTab);
            app.EqualizedImageLabel.HorizontalAlignment = 'center';
            app.EqualizedImageLabel.FontName = 'Yu Gothic Medium';
            app.EqualizedImageLabel.FontSize = 20;
            app.EqualizedImageLabel.Position = [537 430 428 35];
            app.EqualizedImageLabel.Text = 'Equalized Image';

            % Create OriginalImageLabel
            app.OriginalImageLabel = uilabel(app.HistogramTab);
            app.OriginalImageLabel.HorizontalAlignment = 'center';
            app.OriginalImageLabel.FontName = 'Yu Gothic Medium';
            app.OriginalImageLabel.FontSize = 20;
            app.OriginalImageLabel.Position = [36 430 428 35];
            app.OriginalImageLabel.Text = 'Original Image';

            % Create UnequalizedDisplay
            app.UnequalizedDisplay = uiaxes(app.HistogramTab);
            zlabel(app.UnequalizedDisplay, 'Z')
            app.UnequalizedDisplay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.UnequalizedDisplay.XColor = [1 1 1];
            app.UnequalizedDisplay.XTick = [];
            app.UnequalizedDisplay.YColor = [1 1 1];
            app.UnequalizedDisplay.YTick = [];
            app.UnequalizedDisplay.Position = [82 254 333 174];

            % Create EqualizedDisplay
            app.EqualizedDisplay = uiaxes(app.HistogramTab);
            zlabel(app.EqualizedDisplay, 'Z')
            app.EqualizedDisplay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.EqualizedDisplay.XColor = [1 1 1];
            app.EqualizedDisplay.XTick = [];
            app.EqualizedDisplay.YColor = [1 1 1];
            app.EqualizedDisplay.YTick = [];
            app.EqualizedDisplay.Position = [579 254 331 174];

            % Create UnequalizedHistogram
            app.UnequalizedHistogram = uiaxes(app.HistogramTab);
            zlabel(app.UnequalizedHistogram, 'Z')
            app.UnequalizedHistogram.PlotBoxAspectRatio = [2.11855670103093 1 1];
            app.UnequalizedHistogram.FontName = 'Yu Gothic Light';
            app.UnequalizedHistogram.Position = [90 36 320 184];

            % Create EqualizedHistogram
            app.EqualizedHistogram = uiaxes(app.HistogramTab);
            zlabel(app.EqualizedHistogram, 'Z')
            app.EqualizedHistogram.PlotBoxAspectRatio = [2.11340206185567 1 1];
            app.EqualizedHistogram.FontName = 'Yu Gothic Light';
            app.EqualizedHistogram.Position = [587 36 320 184];

            % Create SpatialFilteringTab
            app.SpatialFilteringTab = uitab(app.TabGroup);
            app.SpatialFilteringTab.Title = 'Spatial Filtering';
            app.SpatialFilteringTab.BackgroundColor = [1 0.9098 0.9608];

            % Create SpatialFilterBrowseButton
            app.SpatialFilterBrowseButton = uibutton(app.SpatialFilteringTab, 'push');
            app.SpatialFilterBrowseButton.ButtonPushedFcn = createCallbackFcn(app, @SpatialFilterBrowseButtonPushed, true);
            app.SpatialFilterBrowseButton.FontName = 'Yu Gothic Medium';
            app.SpatialFilterBrowseButton.FontSize = 16;
            app.SpatialFilterBrowseButton.Position = [44 476 89 34];
            app.SpatialFilterBrowseButton.Text = 'Browse';

            % Create BlurredNoisyImageLabel
            app.BlurredNoisyImageLabel = uilabel(app.SpatialFilteringTab);
            app.BlurredNoisyImageLabel.HorizontalAlignment = 'center';
            app.BlurredNoisyImageLabel.FontName = 'Yu Gothic Medium';
            app.BlurredNoisyImageLabel.FontSize = 20;
            app.BlurredNoisyImageLabel.Position = [552 450 428 35];
            app.BlurredNoisyImageLabel.Text = '';

            % Create KernelsizeLabel
            app.KernelsizeLabel = uilabel(app.SpatialFilteringTab);
            app.KernelsizeLabel.FontName = 'Yu Gothic';
            app.KernelsizeLabel.FontSize = 16;
            app.KernelsizeLabel.Position = [93 171 138 28];
            app.KernelsizeLabel.Text = 'Kernel size';

            % Create FactorKLabel
            app.FactorKLabel = uilabel(app.SpatialFilteringTab);
            app.FactorKLabel.FontName = 'Yu Gothic';
            app.FactorKLabel.FontSize = 16;
            app.FactorKLabel.Position = [93 115 87 28];
            app.FactorKLabel.Text = 'Factor K';

            % Create KernelSizeSpinner
            app.KernelSizeSpinner = uispinner(app.SpatialFilteringTab);
            app.KernelSizeSpinner.Step = 2;
            app.KernelSizeSpinner.Limits = [1 99];
            app.KernelSizeSpinner.FontName = 'Yu Gothic';
            app.KernelSizeSpinner.FontSize = 14;
            app.KernelSizeSpinner.Position = [303 171 101 28];
            app.KernelSizeSpinner.Value = 1;

            % Create EnhanceButton
            app.EnhanceButton = uibutton(app.SpatialFilteringTab, 'push');
            app.EnhanceButton.ButtonPushedFcn = createCallbackFcn(app, @EnhanceButtonPushed, true);
            app.EnhanceButton.BackgroundColor = [0.8196 0.8392 0.851];
            app.EnhanceButton.FontName = 'Yu Gothic';
            app.EnhanceButton.FontSize = 14;
            app.EnhanceButton.Visible = 'off';
            app.EnhanceButton.Position = [93 48 118 31];
            app.EnhanceButton.Text = 'Enhance';

            % Create RemoveNoiseButton
            app.RemoveNoiseButton = uibutton(app.SpatialFilteringTab, 'push');
            app.RemoveNoiseButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveNoiseButtonPushed, true);
            app.RemoveNoiseButton.BackgroundColor = [0.8196 0.8392 0.851];
            app.RemoveNoiseButton.FontName = 'Yu Gothic';
            app.RemoveNoiseButton.FontSize = 14;
            app.RemoveNoiseButton.Visible = 'off';
            app.RemoveNoiseButton.Position = [286 47 118 31];
            app.RemoveNoiseButton.Text = 'Remove Noise';

            % Create EnhancedImageLabel
            app.EnhancedImageLabel = uilabel(app.SpatialFilteringTab);
            app.EnhancedImageLabel.HorizontalAlignment = 'center';
            app.EnhancedImageLabel.FontName = 'Yu Gothic Medium';
            app.EnhancedImageLabel.FontSize = 20;
            app.EnhancedImageLabel.Position = [551 219 428 35];
            app.EnhancedImageLabel.Text = '';

            % Create FactorKEditField
            app.FactorKEditField = uieditfield(app.SpatialFilteringTab, 'numeric');
            app.FactorKEditField.RoundFractionalValues = 'on';
            app.FactorKEditField.FontName = 'Yu Gothic';
            app.FactorKEditField.FontSize = 14;
            app.FactorKEditField.Position = [303 119 101 25];

            % Create UnfilteredImageDisplay
            app.UnfilteredImageDisplay = uiaxes(app.SpatialFilteringTab);
            zlabel(app.UnfilteredImageDisplay, 'Z')
            app.UnfilteredImageDisplay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.UnfilteredImageDisplay.XColor = [1 1 1];
            app.UnfilteredImageDisplay.XTick = [];
            app.UnfilteredImageDisplay.YColor = [1 1 1];
            app.UnfilteredImageDisplay.YTick = [];
            app.UnfilteredImageDisplay.Position = [74 258 352 194];

            % Create Blurred_NoisyDisplay
            app.Blurred_NoisyDisplay = uiaxes(app.SpatialFilteringTab);
            zlabel(app.Blurred_NoisyDisplay, 'Z')
            app.Blurred_NoisyDisplay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.Blurred_NoisyDisplay.XColor = [1 1 1];
            app.Blurred_NoisyDisplay.XTick = [];
            app.Blurred_NoisyDisplay.YColor = [1 1 1];
            app.Blurred_NoisyDisplay.YTick = [];
            app.Blurred_NoisyDisplay.Position = [583 259 352 194];

            % Create EnhancedImageDisplay
            app.EnhancedImageDisplay = uiaxes(app.SpatialFilteringTab);
            zlabel(app.EnhancedImageDisplay, 'Z')
            app.EnhancedImageDisplay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.EnhancedImageDisplay.XColor = [1 1 1];
            app.EnhancedImageDisplay.XTick = [];
            app.EnhancedImageDisplay.YColor = [1 1 1];
            app.EnhancedImageDisplay.YTick = [];
            app.EnhancedImageDisplay.Position = [585 23 352 194];

            % Create FourierITab
            app.FourierITab = uitab(app.TabGroup);
            app.FourierITab.Title = 'Fourier I';
            app.FourierITab.BackgroundColor = [1 0.8706 0.6588];

            % Create FourierBrowseButton
            app.FourierBrowseButton = uibutton(app.FourierITab, 'push');
            app.FourierBrowseButton.ButtonPushedFcn = createCallbackFcn(app, @FourierBrowseButtonPushed, true);
            app.FourierBrowseButton.FontName = 'Yu Gothic Medium';
            app.FourierBrowseButton.FontSize = 16;
            app.FourierBrowseButton.Position = [44 476 89 34];
            app.FourierBrowseButton.Text = 'Browse';

            % Create LogDropDown
            app.LogDropDown = uidropdown(app.FourierITab);
            app.LogDropDown.Items = {'With Log', 'Without Log'};
            app.LogDropDown.ValueChangedFcn = createCallbackFcn(app, @LogDropDownValueChanged, true);
            app.LogDropDown.FontName = 'Yu Gothic';
            app.LogDropDown.FontSize = 14;
            app.LogDropDown.BackgroundColor = [0.8196 0.8392 0.851];
            app.LogDropDown.Position = [231 477 135 32];
            app.LogDropDown.Value = 'With Log';

            % Create MagnitudeLabel
            app.MagnitudeLabel = uilabel(app.FourierITab);
            app.MagnitudeLabel.HorizontalAlignment = 'center';
            app.MagnitudeLabel.FontName = 'Yu Gothic Medium';
            app.MagnitudeLabel.FontSize = 20;
            app.MagnitudeLabel.Position = [537 452 428 35];
            app.MagnitudeLabel.Text = 'Magnitude';

            % Create PhaseLabel
            app.PhaseLabel = uilabel(app.FourierITab);
            app.PhaseLabel.HorizontalAlignment = 'center';
            app.PhaseLabel.FontName = 'Yu Gothic Medium';
            app.PhaseLabel.FontSize = 20;
            app.PhaseLabel.Position = [549 207 428 35];
            app.PhaseLabel.Text = 'Phase';

            % Create OriginalImageLabel_2
            app.OriginalImageLabel_2 = uilabel(app.FourierITab);
            app.OriginalImageLabel_2.HorizontalAlignment = 'center';
            app.OriginalImageLabel_2.FontName = 'Yu Gothic Medium';
            app.OriginalImageLabel_2.FontSize = 20;
            app.OriginalImageLabel_2.Position = [39 391 428 35];
            app.OriginalImageLabel_2.Text = 'Original Image';

            % Create OriginalImageDisplay
            app.OriginalImageDisplay = uiaxes(app.FourierITab);
            zlabel(app.OriginalImageDisplay, 'Z')
            app.OriginalImageDisplay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.OriginalImageDisplay.XColor = [1 1 1];
            app.OriginalImageDisplay.XTick = [];
            app.OriginalImageDisplay.YColor = [1 1 1];
            app.OriginalImageDisplay.YTick = [];
            app.OriginalImageDisplay.Position = [43 86 420 287];

            % Create MagnitudeDisplay
            app.MagnitudeDisplay = uiaxes(app.FourierITab);
            zlabel(app.MagnitudeDisplay, 'Z')
            app.MagnitudeDisplay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.MagnitudeDisplay.XColor = [1 1 1];
            app.MagnitudeDisplay.XTick = [];
            app.MagnitudeDisplay.YColor = [1 1 1];
            app.MagnitudeDisplay.YTick = [];
            app.MagnitudeDisplay.Position = [575 255 352 200];

            % Create PhaseDisplay
            app.PhaseDisplay = uiaxes(app.FourierITab);
            zlabel(app.PhaseDisplay, 'Z')
            app.PhaseDisplay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.PhaseDisplay.XColor = [1 1 1];
            app.PhaseDisplay.XTick = [];
            app.PhaseDisplay.YColor = [1 1 1];
            app.PhaseDisplay.YTick = [];
            app.PhaseDisplay.Position = [575 6 352 203];

            % Create FourierIITab
            app.FourierIITab = uitab(app.TabGroup);
            app.FourierIITab.Title = 'Fourier II';
            app.FourierIITab.BackgroundColor = [1 0.851 0.851];

            % Create FourierIIBrowseButton
            app.FourierIIBrowseButton = uibutton(app.FourierIITab, 'push');
            app.FourierIIBrowseButton.ButtonPushedFcn = createCallbackFcn(app, @FourierIIBrowseButtonPushed, true);
            app.FourierIIBrowseButton.FontName = 'Yu Gothic Medium';
            app.FourierIIBrowseButton.FontSize = 16;
            app.FourierIIBrowseButton.Position = [44 476 89 34];
            app.FourierIIBrowseButton.Text = 'Browse';

            % Create FourierKernelsizeLabel
            app.FourierKernelsizeLabel = uilabel(app.FourierIITab);
            app.FourierKernelsizeLabel.FontName = 'Yu Gothic';
            app.FourierKernelsizeLabel.FontSize = 16;
            app.FourierKernelsizeLabel.Position = [93 171 138 28];
            app.FourierKernelsizeLabel.Text = 'Kernel size';

            % Create FourierKernelSizeSpinner
            app.FourierKernelSizeSpinner = uispinner(app.FourierIITab);
            app.FourierKernelSizeSpinner.Step = 2;
            app.FourierKernelSizeSpinner.Limits = [1 99];
            app.FourierKernelSizeSpinner.FontName = 'Yu Gothic';
            app.FourierKernelSizeSpinner.FontSize = 14;
            app.FourierKernelSizeSpinner.Position = [303 171 101 28];
            app.FourierKernelSizeSpinner.Value = 1;

            % Create FilterButton
            app.FilterButton = uibutton(app.FourierIITab, 'push');
            app.FilterButton.ButtonPushedFcn = createCallbackFcn(app, @FilterButtonPushed, true);
            app.FilterButton.BackgroundColor = [0.8196 0.8392 0.851];
            app.FilterButton.FontName = 'Yu Gothic';
            app.FilterButton.FontSize = 14;
            app.FilterButton.Position = [166 95 167 30];
            app.FilterButton.Text = 'Filter';

            % Create FilteredImageLabel
            app.FilteredImageLabel = uilabel(app.FourierIITab);
            app.FilteredImageLabel.HorizontalAlignment = 'center';
            app.FilteredImageLabel.FontName = 'Yu Gothic Medium';
            app.FilteredImageLabel.FontSize = 20;
            app.FilteredImageLabel.Position = [537 452 428 35];
            app.FilteredImageLabel.Text = 'Filtered Image';

            % Create DifferenceLabel
            app.DifferenceLabel = uilabel(app.FourierIITab);
            app.DifferenceLabel.HorizontalAlignment = 'center';
            app.DifferenceLabel.FontName = 'Yu Gothic Medium';
            app.DifferenceLabel.FontSize = 20;
            app.DifferenceLabel.Position = [540 224 428 35];
            app.DifferenceLabel.Text = 'Difference';

            % Create ImageDisplaay
            app.ImageDisplaay = uiaxes(app.FourierIITab);
            zlabel(app.ImageDisplaay, 'Z')
            app.ImageDisplaay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.ImageDisplaay.XColor = [1 1 1];
            app.ImageDisplaay.XTick = [];
            app.ImageDisplaay.YColor = [1 1 1];
            app.ImageDisplaay.YTick = [];
            app.ImageDisplaay.Position = [74 258 352 194];

            % Create FourierFilterDisplay
            app.FourierFilterDisplay = uiaxes(app.FourierIITab);
            zlabel(app.FourierFilterDisplay, 'Z')
            app.FourierFilterDisplay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.FourierFilterDisplay.XColor = [1 1 1];
            app.FourierFilterDisplay.XTick = [];
            app.FourierFilterDisplay.YColor = [1 1 1];
            app.FourierFilterDisplay.YTick = [];
            app.FourierFilterDisplay.Position = [583 259 352 194];

            % Create DifferenceDisplay
            app.DifferenceDisplay = uiaxes(app.FourierIITab);
            zlabel(app.DifferenceDisplay, 'Z')
            app.DifferenceDisplay.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.DifferenceDisplay.XColor = [1 1 1];
            app.DifferenceDisplay.XTick = [];
            app.DifferenceDisplay.YColor = [1 1 1];
            app.DifferenceDisplay.YTick = [];
            app.DifferenceDisplay.Position = [585 23 352 194];

            % Create RemovePatternTab
            app.RemovePatternTab = uitab(app.TabGroup);
            app.RemovePatternTab.Title = 'Remove Pattern';

            % Create RemovePatternButton
            app.RemovePatternButton = uibutton(app.RemovePatternTab, 'push');
            app.RemovePatternButton.ButtonPushedFcn = createCallbackFcn(app, @RemovePatternButtonPushed, true);
            app.RemovePatternButton.FontName = 'Yu Gothic Medium';
            app.RemovePatternButton.FontSize = 16;
            app.RemovePatternButton.Position = [44 476 136 34];
            app.RemovePatternButton.Text = 'Remove Pattern';

            % Create PatternedImageLabel
            app.PatternedImageLabel = uilabel(app.RemovePatternTab);
            app.PatternedImageLabel.HorizontalAlignment = 'center';
            app.PatternedImageLabel.FontName = 'Yu Gothic Medium';
            app.PatternedImageLabel.FontSize = 20;
            app.PatternedImageLabel.Position = [34 393 286 35];
            app.PatternedImageLabel.Text = 'Patterned Image';

            % Create MaskingLabel
            app.MaskingLabel = uilabel(app.RemovePatternTab);
            app.MaskingLabel.HorizontalAlignment = 'center';
            app.MaskingLabel.FontName = 'Yu Gothic Medium';
            app.MaskingLabel.FontSize = 20;
            app.MaskingLabel.Position = [367 393 286 35];
            app.MaskingLabel.Text = 'Masking';

            % Create PatternRemovedLabel
            app.PatternRemovedLabel = uilabel(app.RemovePatternTab);
            app.PatternRemovedLabel.HorizontalAlignment = 'center';
            app.PatternRemovedLabel.FontName = 'Yu Gothic Medium';
            app.PatternRemovedLabel.FontSize = 20;
            app.PatternRemovedLabel.Position = [702 393 286 35];
            app.PatternRemovedLabel.Text = 'Pattern Removed';

            % Create Mask
            app.Mask = uiaxes(app.RemovePatternTab);
            zlabel(app.Mask, 'Z')
            app.Mask.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.Mask.XColor = [1 1 1];
            app.Mask.XTick = [];
            app.Mask.YColor = [1 1 1];
            app.Mask.YTick = [];
            app.Mask.Position = [397 53 240 323];

            % Create ImageWithoutPattern
            app.ImageWithoutPattern = uiaxes(app.RemovePatternTab);
            zlabel(app.ImageWithoutPattern, 'Z')
            app.ImageWithoutPattern.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.ImageWithoutPattern.XColor = [1 1 1];
            app.ImageWithoutPattern.XTick = [];
            app.ImageWithoutPattern.YColor = [1 1 1];
            app.ImageWithoutPattern.YTick = [];
            app.ImageWithoutPattern.Position = [725 53 240 323];

            % Create ImageWithPattern
            app.ImageWithPattern = uiaxes(app.RemovePatternTab);
            zlabel(app.ImageWithPattern, 'Z')
            app.ImageWithPattern.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.ImageWithPattern.XColor = [1 1 1];
            app.ImageWithPattern.XTick = [];
            app.ImageWithPattern.YColor = [1 1 1];
            app.ImageWithPattern.YTick = [];
            app.ImageWithPattern.Position = [57 53 240 323];

            % Create NoiseTab
            app.NoiseTab = uitab(app.TabGroup);
            app.NoiseTab.Title = 'Noise';
            app.NoiseTab.BackgroundColor = [1 0.898 0.8314];

            % Create CreateImageButton
            app.CreateImageButton = uibutton(app.NoiseTab, 'push');
            app.CreateImageButton.ButtonPushedFcn = createCallbackFcn(app, @CreateImageButtonPushed, true);
            app.CreateImageButton.FontName = 'Yu Gothic Medium';
            app.CreateImageButton.FontSize = 16;
            app.CreateImageButton.Position = [43 476 138 34];
            app.CreateImageButton.Text = 'Create Image';

            % Create NoiseTypeDropDown
            app.NoiseTypeDropDown = uidropdown(app.NoiseTab);
            app.NoiseTypeDropDown.Items = {'No Noise', 'Gaussian Noise', 'Uniform Noise'};
            app.NoiseTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @NoiseTypeDropDownValueChanged, true);
            app.NoiseTypeDropDown.FontName = 'Yu Gothic';
            app.NoiseTypeDropDown.FontSize = 14;
            app.NoiseTypeDropDown.BackgroundColor = [0.8196 0.8392 0.851];
            app.NoiseTypeDropDown.Position = [567 374 185 32];
            app.NoiseTypeDropDown.Value = 'No Noise';

            % Create NoiseTypeLabel
            app.NoiseTypeLabel = uilabel(app.NoiseTab);
            app.NoiseTypeLabel.FontName = 'Yu Gothic';
            app.NoiseTypeLabel.FontSize = 16;
            app.NoiseTypeLabel.Position = [403 377 138 28];
            app.NoiseTypeLabel.Text = 'Noise Type';

            % Create ConfirmRegionofInterestButton
            app.ConfirmRegionofInterestButton = uibutton(app.NoiseTab, 'push');
            app.ConfirmRegionofInterestButton.ButtonPushedFcn = createCallbackFcn(app, @ConfirmRegionofInterestButtonPushed, true);
            app.ConfirmRegionofInterestButton.HorizontalAlignment = 'left';
            app.ConfirmRegionofInterestButton.BackgroundColor = [0.8196 0.8392 0.851];
            app.ConfirmRegionofInterestButton.FontName = 'Yu Gothic';
            app.ConfirmRegionofInterestButton.FontSize = 14;
            app.ConfirmRegionofInterestButton.Position = [789 373 185 34];
            app.ConfirmRegionofInterestButton.Text = 'Confirm Region of Interest';

            % Create MeanLabel
            app.MeanLabel = uilabel(app.NoiseTab);
            app.MeanLabel.FontName = 'Yu Gothic';
            app.MeanLabel.FontSize = 16;
            app.MeanLabel.Position = [403 310 120 28];
            app.MeanLabel.Text = 'Mean';

            % Create MeanEditField
            app.MeanEditField = uieditfield(app.NoiseTab, 'text');
            app.MeanEditField.FontName = 'Yu Gothic';
            app.MeanEditField.FontSize = 14;
            app.MeanEditField.Position = [503 311 123 27];

            % Create VarianceLabel
            app.VarianceLabel = uilabel(app.NoiseTab);
            app.VarianceLabel.FontName = 'Yu Gothic';
            app.VarianceLabel.FontSize = 16;
            app.VarianceLabel.Position = [753 307 120 28];
            app.VarianceLabel.Text = 'Variance';

            % Create VarianceEditField
            app.VarianceEditField = uieditfield(app.NoiseTab, 'text');
            app.VarianceEditField.FontName = 'Yu Gothic';
            app.VarianceEditField.FontSize = 14;
            app.VarianceEditField.Position = [853 308 123 27];

            % Create CreatedImage
            app.CreatedImage = uiaxes(app.NoiseTab);
            zlabel(app.CreatedImage, 'Z')
            app.CreatedImage.PlotBoxAspectRatio = [1.06441717791411 1 1];
            app.CreatedImage.YDir = 'reverse';
            app.CreatedImage.XAxisLocation = 'top';
            app.CreatedImage.XColor = [1 1 1];
            app.CreatedImage.XTick = [];
            app.CreatedImage.YColor = [1 1 1];
            app.CreatedImage.YTick = [];
            app.CreatedImage.Position = [35 85 349 357];

            % Create ROIHistogram
            app.ROIHistogram = uiaxes(app.NoiseTab);
            zlabel(app.ROIHistogram, 'Z')
            app.ROIHistogram.PlotBoxAspectRatio = [2.11855670103093 1 1];
            app.ROIHistogram.FontName = 'Yu Gothic Light';
            app.ROIHistogram.Position = [397 85 573 184];

            % Create BackProjectionTab
            app.BackProjectionTab = uitab(app.TabGroup);
            app.BackProjectionTab.Title = 'Back Projection';
            app.BackProjectionTab.BackgroundColor = [0.9843 0.902 1];

            % Create StartButton
            app.StartButton = uibutton(app.BackProjectionTab, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.FontName = 'Yu Gothic Medium';
            app.StartButton.FontSize = 16;
            app.StartButton.Position = [44 476 67 34];
            app.StartButton.Text = 'Start';

            % Create ScheppLoganPhantomLabel
            app.ScheppLoganPhantomLabel = uilabel(app.BackProjectionTab);
            app.ScheppLoganPhantomLabel.HorizontalAlignment = 'center';
            app.ScheppLoganPhantomLabel.FontName = 'Yu Gothic Medium';
            app.ScheppLoganPhantomLabel.FontSize = 20;
            app.ScheppLoganPhantomLabel.Position = [34 393 286 35];
            app.ScheppLoganPhantomLabel.Text = 'Schepp-Logan Phantom';

            % Create SinogramLabel
            app.SinogramLabel = uilabel(app.BackProjectionTab);
            app.SinogramLabel.HorizontalAlignment = 'center';
            app.SinogramLabel.FontName = 'Yu Gothic Medium';
            app.SinogramLabel.FontSize = 20;
            app.SinogramLabel.Position = [367 393 286 35];
            app.SinogramLabel.Text = 'Sinogram';

            % Create LaminogramLabel
            app.LaminogramLabel = uilabel(app.BackProjectionTab);
            app.LaminogramLabel.HorizontalAlignment = 'center';
            app.LaminogramLabel.FontName = 'Yu Gothic Medium';
            app.LaminogramLabel.FontSize = 20;
            app.LaminogramLabel.Position = [702 393 286 35];
            app.LaminogramLabel.Text = 'Laminogram';

            % Create LaminogramDropDown
            app.LaminogramDropDown = uidropdown(app.BackProjectionTab);
            app.LaminogramDropDown.Items = {'Unfiltered Laminogram (0,20,40,60,160)', 'UnFiltered Laminogram (0:179)', 'RamLak', 'Hamming'};
            app.LaminogramDropDown.ValueChangedFcn = createCallbackFcn(app, @LaminogramDropDownValueChanged, true);
            app.LaminogramDropDown.FontName = 'Yu Gothic';
            app.LaminogramDropDown.FontSize = 14;
            app.LaminogramDropDown.BackgroundColor = [0.8196 0.8392 0.851];
            app.LaminogramDropDown.Position = [658 477 310 32];
            app.LaminogramDropDown.Value = 'Unfiltered Laminogram (0,20,40,60,160)';

            % Create ManualButton
            app.ManualButton = uibutton(app.BackProjectionTab, 'push');
            app.ManualButton.ButtonPushedFcn = createCallbackFcn(app, @ManualButtonPushed, true);
            app.ManualButton.BackgroundColor = [0.8196 0.8392 0.851];
            app.ManualButton.FontName = 'Yu Gothic';
            app.ManualButton.FontSize = 14;
            app.ManualButton.Position = [533 478 68 30];
            app.ManualButton.Text = 'Manual';

            % Create BuiltinButton
            app.BuiltinButton = uibutton(app.BackProjectionTab, 'push');
            app.BuiltinButton.ButtonPushedFcn = createCallbackFcn(app, @BuiltinButtonPushed, true);
            app.BuiltinButton.BackgroundColor = [0.8196 0.8392 0.851];
            app.BuiltinButton.FontName = 'Yu Gothic';
            app.BuiltinButton.FontSize = 14;
            app.BuiltinButton.Position = [455 478 68 30];
            app.BuiltinButton.Text = 'Built-in';

            % Create ScheppLoganPhantom
            app.ScheppLoganPhantom = uiaxes(app.BackProjectionTab);
            zlabel(app.ScheppLoganPhantom, 'Z')
            app.ScheppLoganPhantom.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.ScheppLoganPhantom.XColor = [1 1 1];
            app.ScheppLoganPhantom.XTick = [];
            app.ScheppLoganPhantom.YColor = [1 1 1];
            app.ScheppLoganPhantom.YTick = [];
            app.ScheppLoganPhantom.Position = [28 55 298 323];

            % Create Sinogram
            app.Sinogram = uiaxes(app.BackProjectionTab);
            zlabel(app.Sinogram, 'Z')
            app.Sinogram.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.Sinogram.XColor = [1 1 1];
            app.Sinogram.XTick = [];
            app.Sinogram.YColor = [1 1 1];
            app.Sinogram.YTick = [];
            app.Sinogram.Position = [361 54 298 323];

            % Create Laminogram
            app.Laminogram = uiaxes(app.BackProjectionTab);
            zlabel(app.Laminogram, 'Z')
            app.Laminogram.PlotBoxAspectRatio = [1.8974358974359 1 1];
            app.Laminogram.XColor = [1 1 1];
            app.Laminogram.XTick = [];
            app.Laminogram.YColor = [1 1 1];
            app.Laminogram.YTick = [];
            app.Laminogram.Position = [696 54 298 323];

            % Create ImageMorphologyTab
            app.ImageMorphologyTab = uitab(app.TabGroup);
            app.ImageMorphologyTab.Title = 'Image Morphology';
            app.ImageMorphologyTab.BackgroundColor = [0.9216 0.9882 0.9451];

            % Create DisplayButton
            app.DisplayButton = uibutton(app.ImageMorphologyTab, 'push');
            app.DisplayButton.ButtonPushedFcn = createCallbackFcn(app, @DisplayButtonPushed, true);
            app.DisplayButton.FontName = 'Yu Gothic Medium';
            app.DisplayButton.FontSize = 16;
            app.DisplayButton.Position = [44 476 100 34];
            app.DisplayButton.Text = 'Display';

            % Create FingerprintImageLabel
            app.FingerprintImageLabel = uilabel(app.ImageMorphologyTab);
            app.FingerprintImageLabel.HorizontalAlignment = 'center';
            app.FingerprintImageLabel.FontName = 'Yu Gothic Medium';
            app.FingerprintImageLabel.FontSize = 20;
            app.FingerprintImageLabel.Position = [44 430 428 35];
            app.FingerprintImageLabel.Text = 'Fingerprint Image';

            % Create MorphologicalProcessingLabel
            app.MorphologicalProcessingLabel = uilabel(app.ImageMorphologyTab);
            app.MorphologicalProcessingLabel.HorizontalAlignment = 'center';
            app.MorphologicalProcessingLabel.FontName = 'Yu Gothic Medium';
            app.MorphologicalProcessingLabel.FontSize = 20;
            app.MorphologicalProcessingLabel.Position = [533 430 428 35];
            app.MorphologicalProcessingLabel.Text = 'Morphological Processing';

            % Create DilationButton
            app.DilationButton = uibutton(app.ImageMorphologyTab, 'push');
            app.DilationButton.ButtonPushedFcn = createCallbackFcn(app, @DilationButtonPushed, true);
            app.DilationButton.BackgroundColor = [0.8196 0.8392 0.851];
            app.DilationButton.FontName = 'Yu Gothic';
            app.DilationButton.FontSize = 14;
            app.DilationButton.Position = [712 78 69 30];
            app.DilationButton.Text = 'Dilation';

            % Create ErosionButton
            app.ErosionButton = uibutton(app.ImageMorphologyTab, 'push');
            app.ErosionButton.ButtonPushedFcn = createCallbackFcn(app, @ErosionButtonPushed, true);
            app.ErosionButton.BackgroundColor = [0.8196 0.8392 0.851];
            app.ErosionButton.FontName = 'Yu Gothic';
            app.ErosionButton.FontSize = 14;
            app.ErosionButton.Position = [627 78 69 30];
            app.ErosionButton.Text = 'Erosion';

            % Create ClosingButton
            app.ClosingButton = uibutton(app.ImageMorphologyTab, 'push');
            app.ClosingButton.ButtonPushedFcn = createCallbackFcn(app, @ClosingButtonPushed, true);
            app.ClosingButton.BackgroundColor = [0.8196 0.8392 0.851];
            app.ClosingButton.FontName = 'Yu Gothic';
            app.ClosingButton.FontSize = 14;
            app.ClosingButton.Position = [878 78 69 30];
            app.ClosingButton.Text = 'Closing';

            % Create OpeningButton
            app.OpeningButton = uibutton(app.ImageMorphologyTab, 'push');
            app.OpeningButton.ButtonPushedFcn = createCallbackFcn(app, @OpeningButtonPushed, true);
            app.OpeningButton.BackgroundColor = [0.8196 0.8392 0.851];
            app.OpeningButton.FontName = 'Yu Gothic';
            app.OpeningButton.FontSize = 14;
            app.OpeningButton.Position = [796 78 69 30];
            app.OpeningButton.Text = 'Opening';

            % Create Button
            app.Button = uibutton(app.ImageMorphologyTab, 'push');
            app.Button.ButtonPushedFcn = createCallbackFcn(app, @ButtonPushed, true);
            app.Button.BackgroundColor = [0.8196 0.8392 0.851];
            app.Button.FontName = 'Yu Gothic';
            app.Button.FontSize = 14;
            app.Button.Position = [627 30 69 30];
            app.Button.Text = '1';

            % Create Button_2
            app.Button_2 = uibutton(app.ImageMorphologyTab, 'push');
            app.Button_2.ButtonPushedFcn = createCallbackFcn(app, @Button_2Pushed, true);
            app.Button_2.BackgroundColor = [0.8196 0.8392 0.851];
            app.Button_2.FontName = 'Yu Gothic';
            app.Button_2.FontSize = 14;
            app.Button_2.Position = [712 30 69 30];
            app.Button_2.Text = '2';

            % Create Button_3
            app.Button_3 = uibutton(app.ImageMorphologyTab, 'push');
            app.Button_3.ButtonPushedFcn = createCallbackFcn(app, @Button_3Pushed, true);
            app.Button_3.BackgroundColor = [0.8196 0.8392 0.851];
            app.Button_3.FontName = 'Yu Gothic';
            app.Button_3.FontSize = 14;
            app.Button_3.Position = [796 30 69 30];
            app.Button_3.Text = '3';

            % Create Button_4
            app.Button_4 = uibutton(app.ImageMorphologyTab, 'push');
            app.Button_4.ButtonPushedFcn = createCallbackFcn(app, @Button_4Pushed, true);
            app.Button_4.BackgroundColor = [0.8196 0.8392 0.851];
            app.Button_4.FontName = 'Yu Gothic';
            app.Button_4.FontSize = 14;
            app.Button_4.Position = [878 30 69 30];
            app.Button_4.Text = '4';

            % Create RotateWidthLabel_2
            app.RotateWidthLabel_2 = uilabel(app.ImageMorphologyTab);
            app.RotateWidthLabel_2.FontName = 'Yu Gothic';
            app.RotateWidthLabel_2.FontSize = 16;
            app.RotateWidthLabel_2.Position = [503 79 120 28];
            app.RotateWidthLabel_2.Text = 'Operations';

            % Create RotateHeightLabel_2
            app.RotateHeightLabel_2 = uilabel(app.ImageMorphologyTab);
            app.RotateHeightLabel_2.FontName = 'Yu Gothic';
            app.RotateHeightLabel_2.FontSize = 16;
            app.RotateHeightLabel_2.Position = [502 31 120 28];
            app.RotateHeightLabel_2.Text = 'Noise Removal';

            % Create FingerprintImage
            app.FingerprintImage = uiaxes(app.ImageMorphologyTab);
            zlabel(app.FingerprintImage, 'Z')
            app.FingerprintImage.PlotBoxAspectRatio = [1.44587628865979 1 1];
            app.FingerprintImage.XColor = 'none';
            app.FingerprintImage.XTick = [];
            app.FingerprintImage.YColor = 'none';
            app.FingerprintImage.YTick = [];
            app.FingerprintImage.ZColor = [1 1 1];
            app.FingerprintImage.Position = [66 129 384 283];

            % Create MorphologicalImage
            app.MorphologicalImage = uiaxes(app.ImageMorphologyTab);
            zlabel(app.MorphologicalImage, 'Z')
            app.MorphologicalImage.PlotBoxAspectRatio = [1.44587628865979 1 1];
            app.MorphologicalImage.XColor = 'none';
            app.MorphologicalImage.XTick = [];
            app.MorphologicalImage.YColor = 'none';
            app.MorphologicalImage.YTick = [];
            app.MorphologicalImage.ZColor = [1 1 1];
            app.MorphologicalImage.Position = [571 129 384 283];

            % Show the figure after all components are created
            app.ImageProcessing.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = TaskNine_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.ImageProcessing)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ImageProcessing)
        end
    end
end