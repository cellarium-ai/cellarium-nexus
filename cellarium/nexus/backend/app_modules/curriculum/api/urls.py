from rest_framework.routers import DefaultRouter

from nexus.backend.app_modules.curriculum.api.views import CurriculumViewSet

router = DefaultRouter()
router.register("curriculums", CurriculumViewSet, basename="curriculum")

urlpatterns = router.urls 