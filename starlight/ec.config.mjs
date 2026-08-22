import { defineEcConfig } from '@astrojs/starlight/expressive-code';
import ecConfig from 'starlightdown-starlight/ec-config';

// Code block styling and the R output-cell plugin. Expressive Code requires
// this to be a real module at the project root, so it cannot live in
// .starlightdown/. Spread `ecConfig` and override any key to customise.
export default defineEcConfig({ ...ecConfig });
