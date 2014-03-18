#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lv2.h>
#include "Distortion_BigMuff.h"
#include "ModFilters.h"
#include "OverSample.h"

/**********************************************************************************************************************************************************/

#define PLUGIN_URI "http://portalmod.com/plugins/mod-devel/BigMuffPi"
#define TAMANHO_DO_BUFFER 256
enum {IN, OUT_1, TONE, LEVEL, SUSTAIN, PLUGIN_PORT_COUNT};

/**********************************************************************************************************************************************************/

class Distortion
{
public:
    Distortion() {}
    ~Distortion() {}
    static LV2_Handle instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features);
    static void activate(LV2_Handle instance);
    static void deactivate(LV2_Handle instance);
    static void connect_port(LV2_Handle instance, uint32_t port, void *data);
    static void run(LV2_Handle instance, uint32_t n_samples);
    static void cleanup(LV2_Handle instance);
    static const void* extension_data(const char* uri);
    float *in;
    float *out_1;
    float *tone;
    float *level;
    float *sustain;
    
    double *u;
    double *u2;
    double *y;
    double *y2;
    
    double T;
    double SampleRate;
    
    double h1u_1;
    double h1y_1;
    
    double h2u_1;
    double h2y_1;
    double h2u_2;
    double h2y_2;
    double h2u_3;
    double h2y_3;
    
    double h3u_1;
    double h3y_1;
    double h3u_2;
    double h3y_2;
    double h3u_3;
    double h3y_3;
    
    double u_1;
    double y_1;
    
    double h4u_1;
    double h4y_1;
    double h4u_2;
    double h4y_2;
    double h4u_3;
    double h4y_3;
    
    int cont;
    
    double *Sust;
    int nSust;
    
    double SustainMedia_1;
};

/**********************************************************************************************************************************************************/

static const LV2_Descriptor Descriptor = {
    PLUGIN_URI,
    Distortion::instantiate,
    Distortion::connect_port,
    Distortion::activate,
    Distortion::run,
    Distortion::deactivate,
    Distortion::cleanup,
    Distortion::extension_data
};

/**********************************************************************************************************************************************************/

LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
    if (index == 0) return &Descriptor;
    else return NULL;
}

/**********************************************************************************************************************************************************/

LV2_Handle Distortion::instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features)
{
    //printf("--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**--//**sample: %f\n", samplerate );
    Distortion *plugin = new Distortion();
    
    plugin->cont = 0;
    
    plugin->u = (double*)malloc(2*TAMANHO_DO_BUFFER*sizeof(double));
    plugin->y = (double*)malloc(2*TAMANHO_DO_BUFFER*sizeof(double));
    plugin->u2 = (double*)malloc(8*TAMANHO_DO_BUFFER*sizeof(double));
    plugin->y2 = (double*)malloc(8*TAMANHO_DO_BUFFER*sizeof(double));
    
    plugin->nSust = 50;
    plugin->Sust = (double*)calloc(plugin->nSust,sizeof(double));
    
    plugin->SustainMedia_1 = 0.5;
    
    plugin->T = 1/samplerate;
    plugin->SampleRate = samplerate;
    
    plugin->h1u_1 = 0;
    plugin->h1y_1 = 0;
    
    plugin->h2u_1 = 0;
    plugin->h2y_1 = 0;
    plugin->h2u_2 = 0;
    plugin->h2y_2 = 0;
    plugin->h2u_3 = 0;
    plugin->h2y_3 = 0;
    
    plugin->h3u_1 = 0;
    plugin->h3y_1 = 0;
    plugin->h3u_2 = 0;
    plugin->h3y_2 = 0;
    plugin->h3u_3 = 0;
    plugin->h3y_3 = 0;
    
    plugin->u_1 = 0;
    plugin->y_1 = 0;
    
    plugin->h4u_1 = 0;
    plugin->h4y_1 = 0;
    plugin->h4u_2 = 0;
    plugin->h4y_2 = 0;
    plugin->h4u_3 = 0;
    plugin->h4y_3 = 0;
    
    	
    return (LV2_Handle)plugin;
}

/**********************************************************************************************************************************************************/

void Distortion::activate(LV2_Handle instance)
{
    // TODO: include the activate function code here
}

/**********************************************************************************************************************************************************/

void Distortion::deactivate(LV2_Handle instance)
{
    // TODO: include the deactivate function code here
}

/**********************************************************************************************************************************************************/

void Distortion::connect_port(LV2_Handle instance, uint32_t port, void *data)
{
    Distortion *plugin;
    plugin = (Distortion *) instance;

    switch (port)
    {
        case IN:
            plugin->in = (float*) data;
            break;
        case OUT_1:
            plugin->out_1 = (float*) data;
            break;
        case TONE:
            plugin->tone = (float*) data;
            break;
        case LEVEL:
            plugin->level = (float*) data;
            break;
        case SUSTAIN:
            plugin->sustain = (float*) data;
            break;
    }
}

/**********************************************************************************************************************************************************/

void Distortion::run(LV2_Handle instance, uint32_t n_samples)
{
    Distortion *plugin;
    plugin = (Distortion *) instance;
    
    if( (n_samples > TAMANHO_DO_BUFFER) && (plugin->cont == 0) )
	{
		plugin->u = (double*)realloc(plugin->u, 2*n_samples*sizeof(double));
		plugin->y = (double*)realloc(plugin->y, 2*n_samples*sizeof(double));
		plugin->u2 = (double*)realloc(plugin->u2, 8*n_samples*sizeof(double));
		plugin->y2 = (double*)realloc(plugin->y2, 8*n_samples*sizeof(double));
    
		plugin->cont = 1;
    
		return;
	}
    
    double Tone;
    double Level;
    double Sustain;
    double SustainMedia = 0;
    
    Tone = (float)(*(plugin->tone));
    Level = (float)(*(plugin->level));
    Sustain = (float)(*(plugin->sustain));
    
    for (int i=0; i<plugin->nSust-1; i++)
    {
		plugin->Sust[i] = plugin->Sust[i+1];
	}
    plugin->Sust[plugin->nSust-1] = Sustain;
    
    for (int i=0; i<plugin->nSust; i++)
    {
		SustainMedia = SustainMedia + plugin->Sust[i];
	}
	
	SustainMedia = SustainMedia/plugin->nSust;
    
    double T2;
    double T3;
    uint32_t n2;
    uint32_t n3;
    
    for (uint32_t i=1; i<=n_samples; i++)
    {
		plugin->in[i-1] = 5.6234*plugin->in[i-1]; //15dB
	}
	
	//Over 2x
	
	T2 = 0.5*plugin->T;
    Over2(plugin->in, plugin->u, &plugin->h1u_1, n_samples);
    n2 = 2*n_samples;
    
    /*****************************************************************/
    
	Filter1(plugin->u, plugin->y, n2, T2, &plugin->h1u_1, &plugin->h1y_1 );
	
    /*****************************************************************/
   
    for (uint32_t i=1; i<=n2; i++)
    {
		plugin->u[i-1] = plugin->y[i-1];
	}
	
	/*****************************************************************/
	
	Filter2(plugin->u, plugin->y, n2, T2, &plugin->h2u_1, &plugin->h2y_1, &plugin->h2u_2, &plugin->h2y_2, &plugin->h2u_3, &plugin->h2y_3 );
	
   /*****************************************************************/
       
    for (uint32_t i=1; i<=n2; i++)
    {
		plugin->u[i-1] = plugin->y[i-1]; 
	}
	
	/*****************************************************************/
	
	Filter3(plugin->u, plugin->y, n2, T2, &plugin->h3u_1, &plugin->h3y_1, &plugin->h3u_2, &plugin->h3y_2, SustainMedia, plugin->SustainMedia_1 );
	 
	/*****************************************************************/
			
	//Over 4x
	
	T3 = 0.25*T2;
    Over4_Double(plugin->y, plugin->u2, &plugin->u_1, n2);
    n3 = 4*n2;
    
    /*****************************************************************/

	Clip(plugin->u2, plugin->y2, n3, T3, &plugin->u_1, &plugin->y_1);
    
    /*****************************************************************/
    
    for (uint32_t i=1; i<=n3; i++)
    {
		plugin->u2[i-1] = plugin->y2[i-1]; 
	}
	
	/*****************************************************************/
	
	Filter4(plugin->u2, plugin->y2, n3, T3, &plugin->h4u_1, &plugin->h4y_1, &plugin->h4u_2, &plugin->h4y_2, &plugin->h4u_3, &plugin->h4y_3, Tone, Level);
    
    switch (4)
    {
        case 1:
            Down1(plugin->out_1, plugin->y2, n_samples);
            break;
        case 2:
            Down2(plugin->out_1, plugin->y2, n_samples);
            break;
        case 3:
            Down4(plugin->out_1, plugin->y2, n_samples);
            break;
        case 4:
            Down8(plugin->out_1, plugin->y2, n_samples);
            break;
        default:
            Down1(plugin->out_1, plugin->y2, n_samples);
    }
	
	 for (uint32_t i=1; i<=n_samples; i++)
    {
		plugin->out_1[i-1] = plugin->out_1[i-1]/5.6234; //-15dB
	}
	
	plugin->SustainMedia_1 = SustainMedia;
    
}

/**********************************************************************************************************************************************************/

void Distortion::cleanup(LV2_Handle instance)
{
	Distortion *plugin;
	plugin = (Distortion *) instance;
	
	free(plugin->u);
	free(plugin->y);
	free(plugin->u2);
	free(plugin->y2);
	
	free(plugin->Sust);
	
    delete ((Distortion *) instance);
}

/**********************************************************************************************************************************************************/

const void* Distortion::extension_data(const char* uri)
{
    return NULL;
}
