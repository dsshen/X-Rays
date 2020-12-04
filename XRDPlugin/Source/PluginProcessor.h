#pragma once

#include "../JuceLibraryCode/JuceHeader.h"

//==============================================================================

class WavetableOscillator
{
public:
    WavetableOscillator (const AudioSampleBuffer& wavetableToUse)
    : wavetable (wavetableToUse),
    tableSize (wavetable.getNumSamples() - 1)
    {
        jassert (wavetable.getNumChannels() == 1);
    }
    
    void setFrequency (float frequency, float sampleRate)
    {
        float tableSizeOverSampleRate = (float) tableSize / sampleRate;
        tableDelta = frequency * tableSizeOverSampleRate;
        
        // RESET PHASE AND PRAISE THE LORD
        currentIndex = 0.0f;
    }
    
    forcedinline float getNextSample() noexcept
    {
        auto index0 = (unsigned int) currentIndex;
        auto index1 = index0 + 1;
        auto frac = currentIndex - (float) index0;
        auto* table = wavetable.getReadPointer (0);
        auto value0 = table[index0];
        auto value1 = table[index1];
        
        auto currentSample = value0 + frac * (value1 - value0);
        
        if ((currentIndex += tableDelta) >= tableSize)
            currentIndex -= tableSize;
        
        return currentSample;
    }
    
private:
    const AudioSampleBuffer& wavetable;
    const int tableSize;
    float currentIndex = 0.0f, tableDelta = 0.0f;
};

//==============================================================================
/**
*/
class XrdpluginAudioProcessor  : public AudioProcessor,
public Timer,
AudioProcessorValueTreeState::Listener
{
public:
    
    // default stuff from JUCE
public:
    //==============================================================================
    XrdpluginAudioProcessor();
    ~XrdpluginAudioProcessor();
    
    //==============================================================================
    void createWavetable();
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
    
#ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
#endif
    
    void processBlock (AudioBuffer<float>&, MidiBuffer&) override;
    
    //==============================================================================
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;
    
    //==============================================================================
    const String getName() const override;
    
    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;
    
    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;
    
    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;
    
    
    // START IMPORTING YOUR STUFFZ
    //==============================================================================
    //==============================================================================
    //==============================================================================
    

    
    //==============================================================================
    //==============================================================================
    //==============================================================================
    
    // UPDATE FUNCTIONS
    
    // update XRD data
    void updateXRD(double a, double b, double c, double alpha, double beta, double gamma, double lambda)
    {
        isXRDUpdated = false;
        
        int count = 0; // for array indexing
        double maxIntensity = 0.; // for normalizing intensity array
        
        // calculate intensities and 2-theta vals for all hkl using structure factor equation
        for (int h = -1*maxOrder; h <= maxOrder; h++)
        {
            for (int k = -1*maxOrder; k <= maxOrder; k++)
            {
                for (int l = -1*maxOrder; l <= maxOrder; l++)
                {
                    // do math! (see below for details on xrd math functions)
                    double d = get_d(h, k, l, a, b, c, alpha, beta, gamma);
                    double sin_over_lambda = 1 / (2*d);
                    double fNa = get_fNa(sin_over_lambda);
                    double fCl = get_fCl(sin_over_lambda);
                    double Fsqr = get_Fsqr(h, k, l, fNa, fCl);
                    double twoTheta = get_twoTheta(d);
                    double Lp = get_Lp(twoTheta);
                    double intensity = Fsqr * Lp;
                    if (intensity == INFINITY)
                        intensity = 0;       // accounting for hkl = 0,0,0
                    
                    // round twoTheta and intensity to 4 decimal places
                    double twoTheta_round = round(twoTheta*10000.) / 10000.;
                    double intensity_round = round(intensity*10000.) / 10000.;
                    
                    // store twoTheta_round and intensity_round in corresponding arrays
                    twoThetaVals[count] = twoTheta_round;
                    intensityVals[count] = intensity_round;
                    
                    // update max intensity value and increment array index
                    if (intensity > maxIntensity)
                    {
                        maxIntensity = intensity;
                    }
                    count++;
                }
            }
        }
        
        // normalize intensityVals by setting max val to 1
        for (int i = 0; i < numOscs; i++)
        {
            if (!isnan(intensityVals[i]))
            {
                double temp = intensityVals[i] / maxIntensity;
                intensityVals[i] = temp;
            }
            else
                intensityVals[i] = 0.0;
            
            // scale intensity logarithmically
            targetScaledAmps[i] = pow(intensityVals[i], ampExp);
        }
        
        isXRDUpdated = true;
    }
    
    // update frequencies
    void updateFreqs()
    {
        isFreqUpdated = false;
        updateXRD(a_valTimed, b_valTimed, c_valTimed, alpha_valTimed, beta_valTimed, gamma_valTimed, *lambda);
        
        for (int i = 0; i < numOscs; ++i)
        {
            if (!isnan(twoThetaVals[i]))
            {
                targetFrequency[i] = pow(twoThetaVals[i], skewFactor);
            }
            else
            {
                targetFrequency[i] = 0.0;
            }
        }
        
        isFreqUpdated = true;
    }
    
    // Timer callbacks
    void timerCallback() override
    {
        // update xrd variables and make new wavetables every timerCallback() instance
        bool hasChanged = false;
        if (a_valTimed != *a_val)
        {
            a_valTimed = *a_val;
            hasChanged = true;
        }
        if (b_valTimed != *b_val)
        {
            b_valTimed = *b_val;
            hasChanged = true;
        }
        if (c_valTimed != *c_val)
        {
            c_valTimed = *c_val;
            hasChanged = true;
        }
        if (alpha_valTimed != *alpha_val)
        {
            alpha_valTimed = *alpha_val;
            hasChanged = true;
        }
        if (beta_valTimed != *beta_val)
        {
            beta_valTimed = *beta_val;
            hasChanged = true;
        }
        if (gamma_valTimed != *gamma_val)
        {
            gamma_valTimed = *gamma_val;
            hasChanged = true;
        }
//        if (lambdaTimed != *lambda)
//        {
//            lambdaTimed = *lambda;
//            hasChanged = true;
//        }
        if (hasChanged)
        {
            updateFreqs();
            for (auto i = 0; i < numOscs; ++i)
            {
                // set frequencies of unused bank, resetting phase in the process
                if (switchedToBank2)
                {
                    auto osc = oscbank1.getUnchecked(i);
                    osc->setFrequency(targetFrequency[i], currentSampleRate);
                }
                else
                {
                    auto osc = oscbank2.getUnchecked(i);
                    osc->setFrequency(targetFrequency[i], currentSampleRate);
                }
            }
            crossFadeNeeded = true;
        }
    }
    
    //==============================================================================
    //==============================================================================
    //==============================================================================
    
    // XRD FUNCTIONS
    
    // return lattice spacing d based on lattice parameters, general triclinic system
    double get_d(int h, int k, int l, double a, double b, double c, double alpha, double beta, double gamma)
    {
        // convert angles to radians
        double alpha_r = degreesToRadians(alpha);
        double beta_r = degreesToRadians(beta);
        double gamma_r = degreesToRadians(gamma);
        
        // compute and store sines and cosines of angles
        double sin_alpha = std::sin(alpha_r);
        double sin_beta = std::sin(beta_r);
        double sin_gamma = std::sin(gamma_r);
        double cos_alpha = std::cos(alpha_r);
        double cos_beta = std::cos(beta_r);
        double cos_gamma = std::cos(gamma_r);
        
        // compute numerator terms
        double num1 = sin_alpha*sin_alpha*(h*h)/(a*a);
        double num2 = sin_beta*sin_beta*(k*k)/(b*b);
        double num3 = sin_gamma*sin_gamma*(l*l)/(c*c);
        double num4 = cos_alpha*(2*k*l)/(b*c);
        double num5 = cos_beta*(2*h*l)/(a*c);
        double num6 = cos_gamma*(2*k*k)/(a*b);
        
        // computer denominator
        double denom = 1 - (cos_alpha*cos_alpha) - (cos_beta*cos_beta) - (cos_gamma*cos_gamma) + (2*cos_alpha*cos_beta*cos_gamma);
        
        // compute and return d
        double one_over_d_squared = (num1+num2+num3+num4+num5+num6)/denom;
        double d = 1 / std::sqrt(one_over_d_squared);
        return d;
    }
    
    // return atomic scattering factor of Na+ using 6-degree poly fit of raw data from IUCr database
    // x = sin(theta) / lamba (from Bragg equation)
    double get_fNa(double x)
    {
        // store poly coefficients
        double p1, p2, p3, p4, p5, p6, p7;
        p1 = -6.063;
        p2 = 41.23;
        p3 = -107.9;
        p4 = 132.4;
        p5 = -68.93;
        p6 = 1.054;
        p7 = 10.;
        
        // compute polynomial
        double y = p1*(x*x*x*x*x*x) + p2*(x*x*x*x*x) + p3*(x*x*x*x) + p4*(x*x*x) + p5*(x*x) + p6*x + p7;
        return y;
    }
    
    // return atomic scattering factor of Cl- using 6-degree poly fit of raw data from IUCr database
    // x = sin(theta) / lamba (from Bragg equation)
    double get_fCl(double x)
    {
        // store poly coefficients
        double p1, p2, p3, p4, p5, p6, p7;
        p1 = -7.404;
        p2 = 39.97;
        p3 = -72.37;
        p4 = 38.52;
        p5 = 26.79;
        p6 = -40.1;
        p7 = 18.95;
        
        // compute polynomial
        double y = p1*(x*x*x*x*x*x) + p2*(x*x*x*x*x) + p3*(x*x*x*x) + p4*(x*x*x) + p5*(x*x) + p6*x + p7;
        return y;
    }
    
    // compute structure-factor-squared (F^2) for given hkl (NaCl crystal)
    double get_Fsqr(int h, int k, int l, double fNa, double fCl)
    {
        double twoPi = 2*M_PI;
        
        // compute cosine component of F
        double cos_part = fNa * std::cos(twoPi * (h*0.0 + k*0.0 + l*0.0))
        + fNa * std::cos(twoPi * (h*0.0 + k*0.5 + l*0.5))
        + fNa * std::cos(twoPi * (h*0.5 + k*0.0 + l*0.5))
        + fNa * std::cos(twoPi * (h*0.5 + k*0.5 + l*0.0))
        + fCl * std::cos(twoPi * (h*0.5 + k*0.0 + l*0.0))
        + fCl * std::cos(twoPi * (h*0.0 + k*0.5 + l*0.0))
        + fCl * std::cos(twoPi * (h*0.0 + k*0.0 + l*0.5))
        + fCl * std::cos(twoPi * (h*0.5 + k*0.5 + l*0.5));
        
        // compute sine component of F
        double sin_part = fNa * std::sin(twoPi * (h*0.0 + k*0.0 + l*0.0))
        + fNa * std::sin(twoPi * (h*0.0 + k*0.5 + l*0.5))
        + fNa * std::sin(twoPi * (h*0.5 + k*0.0 + l*0.5))
        + fNa * std::sin(twoPi * (h*0.5 + k*0.5 + l*0.0))
        + fCl * std::sin(twoPi * (h*0.5 + k*0.0 + l*0.0))
        + fCl * std::sin(twoPi * (h*0.0 + k*0.5 + l*0.0))
        + fCl * std::sin(twoPi * (h*0.0 + k*0.0 + l*0.5))
        + fCl * std::sin(twoPi * (h*0.5 + k*0.5 + l*0.5));
        
        double Fsqr = cos_part*cos_part + sin_part*sin_part;
        return Fsqr;
    }
    
    // compute 2-theta in degrees, assuming copper anode (lambda = 1.54 angstroms)
    double get_twoTheta(double d)
    {
        // compute 2-theta using Bragg's Law
        double sin_theta = *lambda / (2*d);
        double theta = std::asin(sin_theta);
        double twoTheta_rad = theta*2;
        double twoTheta_deg = radiansToDegrees(twoTheta_rad);
        return twoTheta_deg;
    }
    
    // compute Lorentz-polarization factor (Lp) for given value of 2-theta (in degrees)
    double get_Lp(double twoTheta_deg)
    {
        // get 2-theta and theta in radians
        double twoTheta_rad = degreesToRadians(twoTheta_deg);
        double theta_rad = twoTheta_rad / 2.0;
        
        // store trig quantities in variables
        double cos_theta = std::cos(theta_rad);
        double cos_twoTheta = std::cos(twoTheta_rad);
        double sin_theta = std::sin(theta_rad);
        
        // compute and return Lp
        double Lp = (1.0 + cos_twoTheta * cos_twoTheta) / (2 * (2* sin_theta * sin_theta * cos_theta));
        return Lp;
    }
    
    //==============================================================================
    //==============================================================================
    //==============================================================================
    
    // MIDI PROCESSING FUNCTIONS
    
    // scale controller knob values
    double scale (double val, double min1, double max1, double min2, double max2)
    {
        double relVal = (val - min1) / (max1 - min1);
        double scaledVal = ((max2 - min2) * relVal) + min2;
        return scaledVal;
    }
    
    //==============================================================================
    //==============================================================================
    //==============================================================================
    
    // PARAMETER CHANGE?
    void parameterChanged (const String &parameterID, float newValue) override
    {
        if (parameterID == "attack")
            envelope.setAttackTime(newValue);
        if (parameterID == "decay")
            envelope.setDecayTime(newValue);
        if (parameterID == "sustain")
            envelope.setSustainLevel(newValue);
        if (parameterID == "release")
            envelope.setReleaseTime(newValue);
    }
    
    //==============================================================================
    //==============================================================================
    //==============================================================================
    
    // DECLARATIONS
    
public:
    // adsr envelope
    stk::ADSR envelope;
    float velocityScale = 0.0f;
    
    // for publicly storing parameter values
    float* a_val = nullptr;
    float* b_val = nullptr;
    float* c_val = nullptr;
    float* alpha_val = nullptr;
    float* beta_val = nullptr;
    float* gamma_val = nullptr;
    float* lambda = nullptr;
    float* attack = nullptr;
    float* decay = nullptr;
    float* sustain = nullptr;
    float* release = nullptr;
    AudioProcessorValueTreeState parameters;
    
private:
    // floats for updating xrd vals each timerCallback()
    float a_valTimed;
    float b_valTimed;
    float c_valTimed;
    float alpha_valTimed;
    float beta_valTimed;
    float gamma_valTimed;
//    float lambdaTimed;
    
    // declare arrays for storing 2-theta and intensity vals
    static const int maxOrder = 5; // maximum order of x-ray diffractions
    static const int numOscs = ((2*maxOrder)+1) * ((2*maxOrder)+1) * ((2*maxOrder)+1);
    double twoThetaVals[numOscs] = {0.0};
    double intensityVals[numOscs] = {0.0};
    
    // scaled intensity arrays and objects
    double ampExp = 0.35;
    double targetScaledAmps[numOscs] = {0.0};
    double currentScaledAmps[numOscs] = {0.0};
    
    // smoothing algorithm objects
    double currentFrequency[numOscs] = {0.0};
    double targetFrequency[numOscs] = {0.0};
    
    const int numBufferCycles = 4;  // how many buffers should the crossfade occur over?
    
    double cf = 0.0;
    double cfThisBuffer = 0.0;
    bool isFreqUpdated = false;
    bool isXRDUpdated = false;
    bool isUpdated = false;
    int bufferCycleCounter[numOscs] = {0};
    bool crossFadeNeeded = false;
    bool switchedToBank2 = false;
    
    // skew factor for converting 2theta to Hz
    double skewFactor = 1.73;
    
    // wavetable synthesis objects
    const unsigned int tableSize = 1 << 7;
    float c_level = 0.005f;
    AudioSampleBuffer sineTable;
    OwnedArray<WavetableOscillator> oscbank1;
    OwnedArray<WavetableOscillator> oscbank2;
    double currentSampleRate;
    
    // timer stuff
    int timeInterval = 50;
    
    // midi stuff
    int x; // for storing note number
    
    //==============================================================================
    //==============================================================================
    //==============================================================================
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (XrdpluginAudioProcessor)
};
