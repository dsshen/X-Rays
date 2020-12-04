#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "PluginProcessor.h"

//==============================================================================
/**
 */
class XrdpluginAudioProcessorEditor  : public AudioProcessorEditor,
public Button::Listener
{
public:
    XrdpluginAudioProcessorEditor (XrdpluginAudioProcessor& p, AudioProcessorValueTreeState& vts)
    : AudioProcessorEditor (&p), valueTreeState (vts), processor (p)
    {
        // set up lattice vector sliders and labels
        addAndMakeVisible(a_slider);
        addAndMakeVisible(a_label);
        addAndMakeVisible(b_slider);
        addAndMakeVisible(b_label);
        addAndMakeVisible(c_slider);
        addAndMakeVisible(c_label);
        a_label.setText("a", dontSendNotification);
        a_label.attachToComponent(&a_slider, true);
        b_label.setText("b", dontSendNotification);
        b_label.attachToComponent(&b_slider, true);
        c_label.setText("c", dontSendNotification);
        c_label.attachToComponent(&c_slider, true);
        
        // set up lattice angle sliders and labels
        addAndMakeVisible(alpha_slider);
        addAndMakeVisible(alpha_label);
        addAndMakeVisible(beta_slider);
        addAndMakeVisible(beta_label);
        addAndMakeVisible(gamma_slider);
        addAndMakeVisible(gamma_label);
        alpha_label.setText("alpha", dontSendNotification);
        alpha_label.attachToComponent(&alpha_slider, true);
        beta_label.setText("beta", dontSendNotification);
        beta_label.attachToComponent(&beta_slider, true);
        gamma_label.setText("gamma", dontSendNotification);
        gamma_label.attachToComponent(&gamma_slider, true);
        
        // set up lambda slider
//        addAndMakeVisible(lambda_slider);
//        addAndMakeVisible(lambda_label);
//        lambda_label.setText("lambda", dontSendNotification);
//        lambda_label.attachToComponent(&lambda_slider, true);
        
        // set up ADSR sliders
        addAndMakeVisible(attack_slider);
        addAndMakeVisible(decay_slider);
        addAndMakeVisible(sustain_slider);
        addAndMakeVisible(release_slider);
        attack_label.setText("attack", dontSendNotification);
        attack_label.attachToComponent(&attack_slider, true);
        decay_label.setText("decay", dontSendNotification);
        decay_label.attachToComponent(&decay_slider, true);
        sustain_label.setText("sustain", dontSendNotification);
        sustain_label.attachToComponent(&sustain_slider, true);
        release_label.setText("release", dontSendNotification);
        release_label.attachToComponent(&release_slider, true);

        // set up reset button
        addAndMakeVisible(resetButton);
        resetButton.setButtonText ("RESET");
        resetButton.addListener(this);
        
        // set up debug (note on/off) button
        addAndMakeVisible(debugKey);
        debugKey.setButtonText ("NOTE ON");
        debugKey.addListener(this);
        
        // rounding?
        a_slider.setNumDecimalPlacesToDisplay(2);
        b_slider.setNumDecimalPlacesToDisplay(2);
        c_slider.setNumDecimalPlacesToDisplay(2);
        alpha_slider.setNumDecimalPlacesToDisplay(0);
        beta_slider.setNumDecimalPlacesToDisplay(0);
        gamma_slider.setNumDecimalPlacesToDisplay(0);
//        lambda_slider.setNumDecimalPlacesToDisplay(2);
        
        // attach all slider attachments
        a_attach.reset (new SliderAttachment (valueTreeState, "a_val", a_slider));
        b_attach.reset (new SliderAttachment (valueTreeState, "b_val", b_slider));
        c_attach.reset (new SliderAttachment (valueTreeState, "c_val", c_slider));
        alpha_attach.reset (new SliderAttachment (valueTreeState, "alpha_val", alpha_slider));
        beta_attach.reset (new SliderAttachment (valueTreeState, "beta_val", beta_slider));
        gamma_attach.reset (new SliderAttachment (valueTreeState, "gamma_val", gamma_slider));
//        lambda_attach.reset (new SliderAttachment (valueTreeState, "lambda", lambda_slider));
        attack_attach.reset (new SliderAttachment (valueTreeState, "attack", attack_slider));
        decay_attach.reset (new SliderAttachment (valueTreeState, "decay", decay_slider));
        sustain_attach.reset (new SliderAttachment (valueTreeState, "sustain", sustain_slider));
        release_attach.reset (new SliderAttachment (valueTreeState, "release", release_slider));
        
        // set window size
        setSize (600, 680);
        
        // set listeners for all parameters
//        valueTreeState.addParameterListener("a_val", this);
//        valueTreeState.addParameterListener("b_val", this);
//        valueTreeState.addParameterListener("c_val", this);
//        valueTreeState.addParameterListener("alpha_val", this);
//        valueTreeState.addParameterListener("beta_val", this);
//        valueTreeState.addParameterListener("gamma_val", this);
//        valueTreeState.addParameterListener("lambda", this);
//        valueTreeState.addParameterListener("attack", this);
//        valueTreeState.addParameterListener("decay", this);
//        valueTreeState.addParameterListener("sustain", this);
//        valueTreeState.addParameterListener("release", this);
    }
    
    ~XrdpluginAudioProcessorEditor()
    {
        // delete slider attachments first
        a_attach = nullptr;
        b_attach = nullptr;
        c_attach = nullptr;
        alpha_attach = nullptr;
        beta_attach = nullptr;
        gamma_attach = nullptr;
//        lambda_attach = nullptr;
        attack_attach = nullptr;
        decay_attach = nullptr;
        sustain_attach = nullptr;
        release_attach = nullptr;
        
        // remove listeners
//        valueTreeState.removeParameterListener("a_val", this);
//        valueTreeState.removeParameterListener("b_val", this);
//        valueTreeState.removeParameterListener("c_val", this);
//        valueTreeState.removeParameterListener("alpha_val", this);
//        valueTreeState.removeParameterListener("beta_val", this);
//        valueTreeState.removeParameterListener("gamma_val", this);
//        valueTreeState.removeParameterListener("lambda", this);
//        valueTreeState.removeParameterListener("attack", this);
//        valueTreeState.removeParameterListener("decay", this);
//        valueTreeState.removeParameterListener("sustain", this);
//        valueTreeState.removeParameterListener("release", this);
        resetButton.removeListener(this);
        debugKey.removeListener(this);
    }
    
    //==============================================================================
    
    void paint (Graphics&) override;
    void resized() override;
    
    typedef AudioProcessorValueTreeState::SliderAttachment SliderAttachment;
    typedef AudioProcessorValueTreeState::ButtonAttachment ButtonAttachment;
    
private:
    // headers for virtual callback functions
    //void sliderValueChanged (Slider* slider) override;
    void buttonClicked (Button* button) override;
    void buttonStateChanged (Button* button) override;
//    void parameterChanged (const String &parameterID, float newValue) override;
    
    AudioProcessorValueTreeState& valueTreeState;
    
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    XrdpluginAudioProcessor& processor;
    
    // Slider declarations
    Slider a_slider;
    Slider b_slider;
    Slider c_slider;
    Slider alpha_slider;
    Slider beta_slider;
    Slider gamma_slider;
//    Slider lambda_slider;
    Slider attack_slider;
    Slider decay_slider;
    Slider sustain_slider;
    Slider release_slider;
    
    // Label declarations (for sliders)
    Label a_label;
    Label b_label;
    Label c_label;
    Label alpha_label;
    Label beta_label;
    Label gamma_label;
//    Label lambda_label;
    Label attack_label;
    Label decay_label;
    Label sustain_label;
    Label release_label;
    
    // Attachment declarations
    std::unique_ptr<SliderAttachment> a_attach;
    std::unique_ptr<SliderAttachment> b_attach;
    std::unique_ptr<SliderAttachment> c_attach;
    std::unique_ptr<SliderAttachment> alpha_attach;
    std::unique_ptr<SliderAttachment> beta_attach;
    std::unique_ptr<SliderAttachment> gamma_attach;
//    std::unique_ptr<SliderAttachment> lambda_attach;
    std::unique_ptr<SliderAttachment> attack_attach;
    std::unique_ptr<SliderAttachment> decay_attach;
    std::unique_ptr<SliderAttachment> sustain_attach;
    std::unique_ptr<SliderAttachment> release_attach;
    
    // Reset (and debug) button
    TextButton resetButton;
    TextButton debugKey;
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (XrdpluginAudioProcessorEditor)
};
