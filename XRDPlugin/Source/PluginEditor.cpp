#include "PluginProcessor.h"
#include "PluginEditor.h"


//==============================================================================
//==============================================================================
// make things pretty
void XrdpluginAudioProcessorEditor::paint (Graphics& g)
{
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));
}

void XrdpluginAudioProcessorEditor::resized()
{
    auto sliderLeft = 120;
    auto area = getLocalBounds();
    a_slider.setBounds (sliderLeft, 50, getWidth() - sliderLeft - 10, 20);
    b_slider.setBounds (sliderLeft, 80, getWidth() - sliderLeft - 10, 20);
    c_slider.setBounds (sliderLeft, 110, getWidth() - sliderLeft - 10, 20);
    alpha_slider.setBounds (sliderLeft, 180, getWidth() - sliderLeft - 10, 20);
    beta_slider.setBounds (sliderLeft, 210, getWidth() - sliderLeft - 10, 20);
    gamma_slider.setBounds (sliderLeft, 240, getWidth() - sliderLeft - 10, 20);
//    lambda_slider.setBounds (sliderLeft, 290, getWidth() - sliderLeft - 10, 20);
    attack_slider.setBounds (sliderLeft, 340, getWidth() - sliderLeft - 10, 20);
    decay_slider.setBounds (sliderLeft, 370, getWidth() - sliderLeft - 10, 20);
    sustain_slider.setBounds (sliderLeft, 400, getWidth() - sliderLeft - 10, 20);
    release_slider.setBounds (sliderLeft, 430, getWidth() - sliderLeft - 10, 20);
    resetButton.setBounds (10, 500, getWidth() - 20, 20);
    debugKey.setBounds (10, 530, getWidth() - 20, 20);
}

//==============================================================================
//==============================================================================
//==============================================================================

// SLIDER and BUTTON functions

// update values
//void XrdpluginAudioProcessorEditor::parameterChanged (const String &parameterID, float newValue)
//{
//    if (parameterID == "attack")
//        processor.envelope.setAttackTime(newValue);
//    if (parameterID == "decay")
//        processor.envelope.setDecayTime(newValue);
//    if (parameterID == "sustain")
//        processor.envelope.setSustainLevel(newValue);
//    if (parameterID == "release")
//        processor.envelope.setReleaseTime(newValue);
//}

// button pressed? (mainly for reset button)
void XrdpluginAudioProcessorEditor::buttonClicked (Button* button)
{
    // reset button pressed?
    if (button == &resetButton)
    {
        // NaCl default lattice parameters
        // a, b, c in angstroms; alpha, beta, gamma in degrees
        *processor.a_val = 5.63;
        *processor.b_val = 5.63;
        *processor.c_val = 5.63;
        *processor.alpha_val = 90.;
        *processor.beta_val = 90.;
        *processor.gamma_val = 90.;
//        *processor.lambda = 1.54;
        a_slider.setValue(5.63);
        b_slider.setValue(5.63);
        c_slider.setValue(5.63);
        alpha_slider.setValue(90.);
        beta_slider.setValue(90.);
        gamma_slider.setValue(90.);
//        lambda_slider.setValue(1.54);
    }
}

// debug button
void XrdpluginAudioProcessorEditor::buttonStateChanged (Button* button)
{
    if (button == &debugKey)
    {
        if (debugKey.isDown())
        {
            processor.velocityScale = 1.0;
            processor.envelope.keyOn();
            //std::cout << "on" << std::endl;
        }
        else
        {
            processor.envelope.keyOff();
            //std::cout << "off" << std::endl;
        }
    }
}

//==============================================================================
//==============================================================================
//==============================================================================
