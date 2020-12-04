#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
XrdpluginAudioProcessor::XrdpluginAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     :  AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", AudioChannelSet::stereo(), true)
                     #endif
                       ),
#endif
        parameters (*this, nullptr)
{
    // initialize xrd variables
    parameters.createAndAddParameter ("a_val", "a", String(), Range<float> (1.0f, 10.0f), 5.63f,
                                      [](float value) { return String (value); }, nullptr);
    parameters.createAndAddParameter ("b_val", "b", String(), Range<float> (1.0f, 10.0f), 5.63f,
                                      [](float value) { return String (value); }, nullptr);
    parameters.createAndAddParameter ("c_val", "c", String(), Range<float> (1.0f, 10.0), 5.63f,
                                      [](float value) { return String (value); }, nullptr);
    parameters.createAndAddParameter ("alpha_val", "alpha", String(), Range<float> (0.0f, 180.0f), 90.0f,
                                      [](float value) { return String (value); }, nullptr);
    parameters.createAndAddParameter ("beta_val", "beta", String(), Range<float> (0.0f, 180.0f), 90.0f,
                                      [](float value) { return String (value); }, nullptr);
    parameters.createAndAddParameter ("gamma_val", "gamma", String(), Range<float> (0.0f, 180.0f), 90.0f,
                                      [](float value) { return String (value); }, nullptr);
    parameters.createAndAddParameter ("lambda", "lambda", String(), Range<float> (0.5f, 2.5f), 1.54f,
                                      [](float value) { return String (value); }, nullptr);
    
    // initialize ADSR variables
    parameters.createAndAddParameter ("attack", "attack", String(), Range<float> (0.001f, 1.0f), 0.002f,
                                      [](float value) { return String (value); }, nullptr);
    parameters.createAndAddParameter ("decay", "decay", String(), Range<float> (0.001f, 1.0f), 1.0f,
                                      [](float value) { return String (value); }, nullptr);
    parameters.createAndAddParameter ("sustain", "sustain", String(), Range<float> (0.0f, 1.0f), 0.5f,
                                      [](float value) { return String (value); }, nullptr);
    parameters.createAndAddParameter ("release", "release", String(), Range<float> (0.001f, 1.0f), 0.01f,
                                      [](float value) { return String (value); }, nullptr);
    
    // initialize value tree
    parameters.state = ValueTree (Identifier ("XRDPluginParams"));
    
    // initialize pointers
    a_val = parameters.getRawParameterValue("a_val");
    b_val = parameters.getRawParameterValue("b_val");
    c_val = parameters.getRawParameterValue("c_val");
    alpha_val = parameters.getRawParameterValue("alpha_val");
    beta_val = parameters.getRawParameterValue("beta_val");
    gamma_val = parameters.getRawParameterValue("gamma_val");
    lambda = parameters.getRawParameterValue("lambda");
    attack = parameters.getRawParameterValue("attack");
    decay = parameters.getRawParameterValue("decay");
    sustain = parameters.getRawParameterValue("sustain");
    release = parameters.getRawParameterValue("release");
    
    // initialize timed xrd variables
    a_valTimed = *a_val;
    b_valTimed = *b_val;
    c_valTimed = *c_val;
    alpha_valTimed = *alpha_val;
    beta_valTimed = *beta_val;
    gamma_valTimed = *gamma_val;
    //lambdaTimed = *lambda;
    
    // add listeners for all parameters
    parameters.addParameterListener ("a_val", this);
    parameters.addParameterListener ("b_val", this);
    parameters.addParameterListener ("c_val", this);
    parameters.addParameterListener ("alpha_val", this);
    parameters.addParameterListener ("beta_val", this);
    parameters.addParameterListener ("gamma_val", this);
    parameters.addParameterListener ("lambda", this);
    parameters.addParameterListener ("attack", this);
    parameters.addParameterListener ("decay", this);
    parameters.addParameterListener ("sustain", this);
    parameters.addParameterListener ("release", this);
    
    // start timer
    startTimer(timeInterval);
}

XrdpluginAudioProcessor::~XrdpluginAudioProcessor()
{
    oscbank1.clear(true);
    oscbank2.clear(true);
    a_val = nullptr;
    b_val = nullptr;
    c_val = nullptr;
    alpha_val = nullptr;
    beta_val = nullptr;
    gamma_val = nullptr;
    lambda = nullptr;
    attack = nullptr;
    decay = nullptr;
    sustain = nullptr;
    release = nullptr;
    
    // remove listeners
    parameters.removeParameterListener("a_val", this);
    parameters.removeParameterListener("b_val", this);
    parameters.removeParameterListener("c_val", this);
    parameters.removeParameterListener("alpha_val", this);
    parameters.removeParameterListener("beta_val", this);
    parameters.removeParameterListener("gamma_val", this);
    parameters.removeParameterListener("lambda", this);
    parameters.removeParameterListener("attack", this);
    parameters.removeParameterListener("decay", this);
    parameters.removeParameterListener("sustain", this);
    parameters.removeParameterListener("release", this);
}

//==============================================================================
const String XrdpluginAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool XrdpluginAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool XrdpluginAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool XrdpluginAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double XrdpluginAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int XrdpluginAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int XrdpluginAudioProcessor::getCurrentProgram()
{
    return 0;
}

void XrdpluginAudioProcessor::setCurrentProgram (int index)
{
}

const String XrdpluginAudioProcessor::getProgramName (int index)
{
    return {};
}

void XrdpluginAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void XrdpluginAudioProcessor::createWavetable()
{
    sineTable.setSize (1, tableSize + 1);
    auto* samples = sineTable.getWritePointer (0); // [3]
    auto angleDelta = MathConstants<double>::twoPi / (double) (tableSize - 1); // [4]
    auto currentAngle = 0.0;
    for (auto i = 0; i < tableSize; ++i)
    {
        auto sample = std::sin (currentAngle);     // [5]
        samples[i] = (float) sample;
        currentAngle += angleDelta;
    }
    samples[tableSize] = samples[0];
}

void XrdpluginAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // update freqs for the first time
    updateFreqs();
    
    // initialize ADSR envelope
    envelope.setAllTimes(*attack, *decay, *sustain, *release);
    
    // create wavetable
    createWavetable();
    
    // Prepare to play
    currentSampleRate = sampleRate;
    for (auto i = 0; i < numOscs; ++i)
    {
        // initialize current and target frequencies
        auto* oscillator_1 = new WavetableOscillator (sineTable);
        auto* oscillator_2 = new WavetableOscillator (sineTable);
        if (!isnan(twoThetaVals[i]))
        {
            targetFrequency[i] = pow(twoThetaVals[i], skewFactor);
            currentFrequency[i] = targetFrequency[i];
        }
        else
        {
            targetFrequency[i] = 0.0;
            currentFrequency[i] = 0.0;
        }
        
        // initialize current and target intensities
        currentScaledAmps[i] = targetScaledAmps[i];
        
        oscillator_1->setFrequency (targetFrequency[i], sampleRate);
        oscillator_2->setFrequency (targetFrequency[i], sampleRate);
        oscbank1.add (oscillator_1);
        oscbank2.add (oscillator_2);
    }
}

void XrdpluginAudioProcessor::releaseResources()
{
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool XrdpluginAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void XrdpluginAudioProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    // In case we have more outputs than inputs, this code clears any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    // This is here to avoid people getting screaming feedback
    // when they first compile a plugin, but obviously you don't need to keep
    // this code if your algorithm always overwrites all the output channels.
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

    // This is the place where you'd normally do the guts of your plugin's
    // audio processing...
    // Make sure to reset the state if your inner loop is processing
    // the samples and the outer loop is handling the channels.
    // Alternatively, you can process the samples with the channels
    // interleaved by keeping the same state.
    
    //==============================================================================
    // AUDIO PROCESSING
    
    for (int channel = 0; channel < totalNumOutputChannels; ++channel)
    {
        auto bufferSize = buffer.getNumSamples();
        auto* channelData = buffer.getWritePointer (channel);
        buffer.clear(channel, 0, bufferSize);
        
        // make array for storing adsr envelope values for current audio block
        double env[bufferSize];
        for (int i = 0; i < bufferSize; ++i)
        {
            env[i] = envelope.lastOut();
            envelope.tick();
        }
        
        // for every osc...
        for (int i = 0; i < numOscs; ++i)
        {
            isUpdated = isFreqUpdated && isXRDUpdated;
            
            // if xrd values have changed...
            if (crossFadeNeeded && isUpdated)
            {
                // set up crossfade incrementer
                double cfInc = 1. / (bufferSize * numBufferCycles);
                
                // make pointers to relevant osc objects
                auto* oscillator1 = oscbank1.getUnchecked(i);
                auto* oscillator2 = oscbank2.getUnchecked(i);
                
                // fill in buffer
                for (auto sample = 0; sample < bufferSize; ++sample)
                {
                    float fadeDn = 0.0;
                    float fadeUp = 0.0;
                    if (switchedToBank2)
                    {
                        fadeDn = oscillator2->getNextSample() * currentScaledAmps[i] * c_level * velocityScale * env[sample] * (1.0 - cf);
                        fadeUp = oscillator1->getNextSample() * targetScaledAmps[i] * c_level * velocityScale * env[sample] * cf;
                    }
                    else
                    {
                        fadeDn = oscillator1->getNextSample() * currentScaledAmps[i] * c_level * velocityScale * env[sample] * (1.0 - cf);
                        fadeUp = oscillator2->getNextSample() * targetScaledAmps[i] * c_level * velocityScale * env[sample] * cf;
                    }
                    
                    channelData[sample] += (fadeDn + fadeUp);
                    
                    cf += cfInc;
                    
                    // at end of buffer...
                    if (sample >= bufferSize - 1)
                    {
                        // if still working on oscs, reset cf
                        if (i < numOscs - 1)
                            cf = cfThisBuffer;
                        
                        // if all oscs covered, update cfThisBuffer
                        else
                            cfThisBuffer = cf;
                    }
                }
                
                ++bufferCycleCounter[i];
                
                // check if we're done crossfading
                if (bufferCycleCounter[i] >= numBufferCycles)
                {
                    bufferCycleCounter[i] = 0;
                    currentScaledAmps[i] = targetScaledAmps[i];
                    if (i >= numOscs - 1)
                    {
                        crossFadeNeeded = false; // stop crossfade when last i has been completed
                        switchedToBank2 = !switchedToBank2; // switch which oscbank is "active"
                        cf = 0.0;
                        cfThisBuffer = 0.0;
                    }
                }
            }
            
            // if you're not moving a knob (i.e. idle)
            else
            {
                for (auto sample = 0; sample < bufferSize; ++sample)
                {
                    float levelSample;
                    if (switchedToBank2)
                    {
                        auto* oscillator2 = oscbank2.getUnchecked(i);
                        levelSample = oscillator2->getNextSample() * currentScaledAmps[i] * c_level * velocityScale * env[sample];
                    }
                    else
                    {
                        auto* oscillator1 = oscbank1.getUnchecked(i);
                        levelSample = oscillator1->getNextSample() * currentScaledAmps[i] * c_level * velocityScale * env[sample];
                    }
                    channelData[sample]  += levelSample;
                }
            }
        }
    }
    
    //==============================================================================
    // MIDI PROCESSING
    
    MidiBuffer processedMidi;
    int time;
    MidiMessage m;
    for (MidiBuffer::Iterator i (midiMessages); i.getNextEvent (m, time);)
    {
        // note on and off
        if (m.isNoteOn())
        {
            int velocity = m.getVelocity();
            velocityScale = velocity / 127.;
            
            // convert midi note to lambda
            x = m.getNoteNumber();
            if ((x > 35) && (x < 73))       // restricted to 4-oct range (upper 3 are on MPKMini)
                                            // additional lower oct for bass hits
            {
                if (x == 60) *lambda = 1.54; // hard-code middle C to copper K-alpha wavelength
                
                else
                {
                    double p1 = 0.000006351;
                    double p2 = -0.0003717;
                    double p3 = 0.02609;
                    double p4 = -0.05929;
                    
                    double lambdaDouble = p1*x*x*x + p2*x*x + p3*x + p4;
                    *lambda = (float) lambdaDouble;
                }
                
                // update freqs and oscs
                updateFreqs();
                for (auto i = 0; i < numOscs; ++i)
                {
                    // set frequencies of ACTIVE bank, resetting phase in the process (no cross-fades here!)
                    if (switchedToBank2)
                    {
                        auto osc = oscbank2.getUnchecked(i);
                        osc->setFrequency(targetFrequency[i], currentSampleRate);
                    }
                    else
                    {
                        auto osc = oscbank1.getUnchecked(i);
                        osc->setFrequency(targetFrequency[i], currentSampleRate);
                    }
                }
            }
            envelope.keyOn();
        }
        if (m.isNoteOff() && (m.getNoteNumber() == x))
        {
            envelope.keyOff();
        }
        processedMidi.addEvent (m, time);
    }
    midiMessages.swapWith (processedMidi);
}

//==============================================================================
bool XrdpluginAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* XrdpluginAudioProcessor::createEditor()
{
    return new XrdpluginAudioProcessorEditor (*this, parameters);
}

//==============================================================================
void XrdpluginAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
    
    auto state = parameters.copyState();
    std::unique_ptr<XmlElement> xml (state.createXml());
    copyXmlToBinary (*xml, destData);
}

void XrdpluginAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
    
    std::unique_ptr<XmlElement> xmlState (getXmlFromBinary (data, sizeInBytes));
    if (xmlState.get() != nullptr)
        if (xmlState->hasTagName (parameters.state.getType()))
        {
            parameters.replaceState (ValueTree::fromXml (*xmlState));
            
            // for some reason ADSR doesn't update w/o this addedline
            envelope.setAllTimes(*attack, *decay, *sustain, *release);
        }
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new XrdpluginAudioProcessor();
}

