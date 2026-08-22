import { defineCollection } from 'astro:content';
import { docsLoader } from '@astrojs/starlight/loaders';
import { docsSchema } from '@astrojs/starlight/schema';
import { sdSchema } from 'starlightdown-starlight/schema';

// `sdSchema` validates the `sd:` frontmatter block that starlightdown writes.
// A malformed field fails the build with the offending file and key named.
export const collections = {
	docs: defineCollection({
		// Astro's default ID generator strips dots, which breaks public S3 method
		// routes such as `predict.ref_fit`. Preserve the source path instead.
		loader: docsLoader({
			generateId: ({ entry }) =>
				entry
					.replace(/\\/g, '/')
					.replace(/\.(?:markdown|mdown|mkdn|mkd|mdwn|md|mdx)$/i, '')
					.replace(/\/index$/, ''),
		}),
		schema: docsSchema({ extend: sdSchema }),
	}),
};
